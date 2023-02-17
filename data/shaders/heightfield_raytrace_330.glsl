#version 330 core

#ifdef VERTEX_SHADER
void main(void)
{
	vec4 vertices[4] = vec4[4](vec4(-1.0, -1.0, 1.0, 1.0),
                               vec4( 1.0, -1.0, 1.0, 1.0),
                               vec4(-1.0,  1.0, 1.0, 1.0),
                               vec4( 1.0,  1.0, 1.0, 1.0));
    vec4 pos = vertices[gl_VertexID];
    gl_Position = pos;
} 
#endif

#ifdef FRAGMENT_SHADER
// Camera data
uniform vec3 CamPos;
uniform vec3 CamLookAt;
uniform vec3 CamUp;
uniform float CamAngleOfViewV;
uniform vec2 iResolution;

// Heightfield data
uniform vec2 a;
uniform vec2 b;
uniform vec2 zRange;
uniform float K;
uniform ivec2 texSize;
uniform int shadingMode;

// Textures: elevation and albedo data
uniform sampler2D heightfield;
uniform sampler2D albedo;

// Raymarching data
uniform int STEPS = 512;
uniform float epsilon = 0.1f;
uniform vec3 lightDir = vec3(-1.0f, -1.0f, -1.3f);

out vec4 Fragment;

// Utility function
float Remap(float x, float oldMin, float oldMax, float newMin, float newMax) {
	return newMin + (newMax - newMin) * ((x - oldMin) / (oldMax - oldMin));
}

// Read height from the heightfield texture given a world point
// p: world position
// returns height at point
float Height(vec2 p) {
	vec2 q = p - a;
	vec2 d = b - a;
	vec2 uv = q / d;
	return Remap(texture(heightfield, uv).r, 0.0f, 1.0f, zRange.x, zRange.y);
}

// Read color from the albedo texture given a world point
// returns color at point
vec4 Albedo(vec2 p) {
	vec2 q = p - a;
	vec2 d = b - a;
	vec2 uv = q / d;
	return texture(albedo, uv);
}

// Compute the 3D gradient using central difference
// p: point
// eps: epsilon for the computation
vec2 Gradient(vec3 p, vec2 eps) {
    vec3 e = vec3(eps.x, eps.y, 0.0);
    return vec2(
        (Height(p.xy + e.xz) - Height(p.xy - e.xz)) / (2.0f * e.x),
        (Height(p.xy + e.zy) - Height(p.xy - e.zy)) / (2.0f * e.y)
    );
}

// Compute the 3D normal using central difference
// p: point
// eps: epsilon for the computation
vec3 Normal(vec3 p, vec2 eps) {
    vec3 e = vec3(eps.x, eps.y, 0.0);
    return normalize(vec3(
        Height(p.xy + e.xz) - Height(p.xy - e.xz),
        Height(p.xy + e.zy) - Height(p.xy - e.zy),
		length(eps.xy)
    ));
}

// Build the ray direction from the fragment position
vec3 BuildRd() {	
	vec3 view = normalize(CamLookAt - CamPos);
	vec3 horizontal = normalize(cross(view, CamUp));
	vec3 vertical = normalize(cross(horizontal, view));
	
	float length = 1.0f;
	float rad = CamAngleOfViewV;
	
	float vLength = tan(rad / 2.0f) * length;
	float hLength = vLength * iResolution.x / iResolution.y;
	
	vertical *= vLength;
	horizontal *= hLength;

	vec2 p = gl_FragCoord.xy;

	p.x = p.x - iResolution.x / 2.0f;
	p.y = p.y - iResolution.y / 2.0f;

	p.x -= 0.5f;
	p.y -= 0.5f;

	p.x /= iResolution.x / 2.0f;
	p.y /= iResolution.y / 2.0f;
	
	return normalize(view * length + horizontal * p.x + vertical * p.y);
}

// Intersection
// va, vb: Signed distance field values
float Intersection(float va, float vb) {
	return max(va, vb);
}

// Box signed distance field 
// p: Point
// va, vb: Vertices of the box
float Box(vec3 p, vec3 va, vec3 vb) {
	vec3 c = 0.5 * (va + vb);
	vec3 r = 0.5 * (vb - va);
	vec3 q = abs(p - c) - r;
	float d = length(max(q, 0.0)) + min(max(q.x, max(q.y, q.z)),0.0);
	return d;
}

// Box signed distance field 
// p: Point
// va, vb: Vertices of the box
float Box(vec2 p, vec2 va, vec2 vb) {
	vec2 c = 0.5 * (va + vb);
	vec2 r = 0.5 * (vb - va);
	vec2 q = abs(p - c) - r;
	float d = length(max(q, 0.0)) + min(max(q.x, q.y), 0.0);
	return d;
}

// Signed distance field object
// p: Point
// returns signed distance value for the scene.
float Map(vec3 p) {
	float t = p.z - Height(p.xy);
	float delta = 0.1f * (zRange.y - zRange.x);
	return Intersection(Box(p, vec3(a.x, a.y, zRange.x - delta), vec3(b.x, b.y, zRange.y + delta)), t);
}

// Intersection test between a ray and the slightly inflated terrain 3D bounding box.
// ro: ray origin
// rd: ray direction
// tN, tF: returned intersection parameterizations
bool IntersectBox(vec3 ro, vec3 rd, out float tN, out float tF) {
	vec3 rinvDir = 1.0 / rd;
	float delta = 0.1 * (zRange.y - zRange.x);
	vec3 tbot = rinvDir * (vec3(a.x, a.y, zRange.x - delta) - ro);
	vec3 ttop = rinvDir * (vec3(b.x, b.y, zRange.y + delta) - ro);
	vec3 tmin = min(ttop, tbot);
	vec3 tmax = max(ttop, tbot);
	vec2 t = max(tmin.xx, tmin.yz);
	float t0 = max(t.x, t.y);
	t = min(tmax.xx, tmax.yz);
	float t1 = min(t.x, t.y);
	tN = t0;
	tF = t1;
	return t1 > max(t0, 0.0);
}

// Sphere tracing of the properly defined signed distance field of the terrain along a ray.
// ro: ray origin
// rd: ray direction
// p: returned intersection point on terrain
// t: returned parameterization along the ray
// s: returned sphere tracing steps
// returns true of intersection occurred, false otherwise.
bool SphereTrace(vec3 ro, vec3 rd, out vec3 p, out float t, out int s) {
    t = 0.0f;
	s = 0;
	
	// Lipschitz bound is dependent on ray direction
	float uz = abs(rd.z);
	float kr = uz + K * sqrt(1.0f - (uz * uz));

	// Check if ray intersects the bounding box of the heightfield
	float ta, tb;
	if (!IntersectBox(ro, rd, ta, tb))
		return false;
	
	// Raymarch
	t = max(ta, 0.0f);
	float d = 0.0f;
    for(s = 0; s < STEPS; s++) 
	{
		if (t > tb)
		    return false;
        p = ro + rd * t;
        d = Map(p);
        if(d <= 0.0) 
			return true;
        t += max(d / kr, epsilon);
    }
    return false;
}

// Compute sky color 
// d: Ray direction
// returns sky color
vec4 ShadeSkyBlue(in vec3 d) {
	vec3 lig = normalize(vec3(0.3, 0.5, 0.6));
	float sun = 0.5 * (dot(lig, d) + 1.0);
	vec3 color = vec3(0.65, 0.85, 1.0) * (0.75 - 0.25 * d.z);
	color += vec3(0.95, 0.90, 0.95) * pow(sun, 18.0);
	return vec4(color, 1.0);
}

// Compute terrain color
// p: intersection point on terrain
// returns terrain color.
vec4 ShadeTerrain(vec3 p) {
	// Terrain sides and bottom
	if (abs(Box(p.xy, a, b)) < epsilon || abs(p.z - zRange.x + 0.1f * (zRange.y - zRange.x)) < epsilon)
		return vec4(0.3f, 0.29f, 0.31f, 1.0f);
	
	// Terrain interior
	if (shadingMode == 0)
	{
		vec3 n = Normal(p, (b - a) / texSize);
		return vec4(0.2 * (vec3(3.0) + 2.0 * n.xyz), 1.0);
	}
	else if (shadingMode == 1)
	{
		vec4 col = Albedo(p.xy);
		vec3 light = -normalize(lightDir);
		vec3 viewDir = normalize(CamPos - p);
		vec3 n = Normal(p, (b - a) / texSize);
		float irradiance = max(dot(light, n), 0.0f) * 1.0f;
		return vec4(col.xyz * irradiance, 1.0f);
	}
	else
		return vec4(1.0, 1.0, 1.0, 1.0);
}

void main() {
	// Compute ray
	vec3 ro = CamPos;
	vec3 rd = BuildRd();

	// Compute Intersection with terrain
	vec4 c;
	vec3 p;
	float t;
	int s;
	bool hit = SphereTrace(ro, rd, p, t, s);
	
	// Shade either sky or terrain
	if (hit) {
		c = ShadeTerrain(p);
	}
	else {
		c = ShadeSkyBlue(rd);
		c = mix(c, vec4(0.85, 0.95, 1.0, 1.0), gl_FragCoord.x * 0.0001 + 0.95 * gl_FragCoord.y / iResolution.y);
	}

	Fragment = c;
}

#endif
