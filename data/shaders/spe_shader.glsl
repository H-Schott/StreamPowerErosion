#version 450 core

#ifdef COMPUTE_SHADER

// In
layout(binding = 0, std430) readonly buffer InElevation { float hf[]; };
layout(binding = 1, std430) readonly buffer InStreamArea { float stream[]; };

// Out
layout(binding = 2, std430) writeonly buffer OutElevation { float out_hf[]; };
layout(binding = 3, std430) writeonly buffer OutStreamArea { float out_stream[]; };

layout(binding = 4, std430) readonly buffer Uplift { float upliftBuffer[]; };



uniform int nx;
uniform int ny;
uniform vec2 a;
uniform vec2 b;
uniform vec2 cellDiag;

// 0: Stream power
// 1: Stream power + Hillslope (Laplacian)
// 2: Stream power + Hillslope (Laplacian) + Debris slope
uniform int erosionMode = 0;

uniform float uplift = 0.01f;
uniform float k = 0.0005f;
uniform float k_d = 10.0f;
uniform float k_h = 2.0f;
uniform float p_sa = 0.8f;
uniform float p_sl = 2.0f;
uniform float dt = 50.0f;

const ivec2 next8[8] = ivec2[8](ivec2(0, 1), ivec2(1, 1), ivec2(1, 0), ivec2(1, -1),
                                ivec2(0, -1), ivec2(-1, -1), ivec2(-1, 0), ivec2(-1, 1));

int ToIndex1D(int i, int j) { return i + nx * j; }

int ToIndex1D(ivec2 p) { return p.x + nx * p.y; }

float Height(ivec2 p) { return hf[ToIndex1D(p)]; }

float UpliftAt(int i, int j) { return upliftBuffer[ToIndex1D(i, j)]; }

vec2 ArrayPoint(ivec2 p) { return a.xy + vec2(p) * cellDiag; }

vec3 Point3D(ivec2 p) {
    return vec3(ArrayPoint(p), Height(p));
}

vec4 Read(ivec2 p) {
    if (p.x < 0 || p.x >= nx || p.y < 0 || p.y >= ny)
        return vec4(0);
    int id = ToIndex1D(p);
    vec4 ret;
    ret.x = hf[id];                // Bedrock elevation
    ret.y = stream[id];            // Stream area
    ret.z = upliftBuffer[id];      // Uplift factor
    return ret;
}

void Write(int id, vec4 data) {
    out_hf[id] = data.x;
    out_stream[id] = data.y;
}

float Slope(ivec2 p, ivec2 q) {
    if (p.x < 0 || p.x >= nx || p.y < 0 || p.y >= ny) return 0.0f;
    if (q.x < 0 || q.x >= nx || q.y < 0 || q.y >= ny) return 0.0f;
    if (p == q) return 0.0f;
    int index_p = ToIndex1D(p.x, p.y);
    int index_q = ToIndex1D(q.x, q.y);
    float d = length(ArrayPoint(q) - ArrayPoint(p));
    return (hf[index_q] - hf[index_p]) / d;
}

float Stream(ivec2 p) {
    if (p.x < 0 || p.x >= nx || p.y < 0 || p.y >= ny) return 0.0f;
    int index_p = ToIndex1D(p.x, p.y);
    return stream[index_p];
}


float Laplacian(ivec2 p) {
    float lapl = 0.0f;
    int i = p.x;
    int j = p.y;

    if (i == 0)
        lapl += (hf[ToIndex1D(i, j)] - 2.0f * hf[ToIndex1D(i + 1, j)] + hf[ToIndex1D(i + 2, j)]) / (cellDiag.x * cellDiag.x);
    else if (i == nx - 1)
        lapl += (hf[ToIndex1D(i, j)] - 2.0f * hf[ToIndex1D(i - 1, j)] + hf[ToIndex1D(i - 2, j)]) / (cellDiag.x * cellDiag.x);
    else
        lapl += (hf[ToIndex1D(i + 1, j)] - 2.0f * hf[ToIndex1D(i, j)] + hf[ToIndex1D(i - 1, j)]) / (cellDiag.x * cellDiag.x);

    if (j == 0)
        lapl += (hf[ToIndex1D(i, j)] - 2.0f * hf[ToIndex1D(i, j + 1)] + hf[ToIndex1D(i, j + 2)]) / (cellDiag.y * cellDiag.y);
    else if (j == ny - 1)
        lapl += (hf[ToIndex1D(i, j)] - 2.0f * hf[ToIndex1D(i, j - 1)] + hf[ToIndex1D(i, j - 2)]) / (cellDiag.y * cellDiag.y);
    else
        lapl += (hf[ToIndex1D(i, j + 1)] - 2.0f * hf[ToIndex1D(i, j)] + hf[ToIndex1D(i, j - 1)]) / (cellDiag.y * cellDiag.y);

    return lapl;
}

ivec2 GetFlowSteepest(ivec2 p) {
    ivec2 d = ivec2(0, 0);
    float maxSlope = 0.0f;
    for (int i = 0; i < 8; i++) {
        float ss = Slope(p + next8[i], p);
        if (ss > maxSlope) {
            maxSlope = ss;
            d = next8[i];
        }
    }
    return d;
}

float WaterSteepest(ivec2 p) {
    float water = 0.0f;
    for (int i = 0; i < 8; i++) {
        ivec2 q = p + next8[i];
        ivec2 fd = GetFlowSteepest(q);
        if (q + fd == p) {
            water += Stream(q);
        }
    }
    return water;
}

layout(local_size_x = 8, local_size_y = 8, local_size_z = 1) in;
void main() {
    int x = int(gl_GlobalInvocationID.x);
    int y = int(gl_GlobalInvocationID.y);
    if (x < 0)   return;
    if (y < 0)   return;
    if (x >= nx) return;
    if (y >= ny) return;

    int id = ToIndex1D(x, y);
    ivec2 p = ivec2(x, y);
    vec4 data = Read(p);

    // Border nodes are fixed to zero (elevation and drainage)
    if (p.x == 0 || p.x == nx - 1 || p.y == 0 || p.y == ny - 1)
    {
        data.x = 0.0f;
        data.y = 1.0f * length(cellDiag);
        Write(id, data);
        return;
    }

    // Flows accumulation at p
    float waterIncr = WaterSteepest(p);

    data.y = 1.0f * length(cellDiag);
    data.y += waterIncr;
   
    // Erosion at p (relative to steepest)
    ivec2 d = GetFlowSteepest(p);
    vec4 receiver = Read(p + d);
    float pslope = abs(Slope(p + d, p));

    float spe = k * pow(data.y, p_sa) * pow(pslope, p_sl);

    float newH = data.x;
    if (erosionMode == 0)       // Stream power
        newH -= dt * (spe);
    else if (erosionMode == 1)  // Stream power + Hillslope erosion (Laplacian)
        newH -= dt * (spe - k_h * Laplacian(p));
    else if (erosionMode == 2)  // Stream power + Hillslope erosion (Laplacian) + Debris flow
        newH -= dt * (spe - k_h * Laplacian(p) - k_d * pslope);
    newH = max(newH, receiver.x);
    newH += dt * uplift * data.z;

    data.x = newH;
    Write(id, data);
}

#endif
