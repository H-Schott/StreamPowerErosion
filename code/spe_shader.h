#ifndef __SPE_SHADER__
#define __SPE_SHADER__

#include <GL/glew.h>
#include <vector>
#include <string>

#include "scalarfield2.h"
#include "shader-api.h"

class GPU_SPE {
private:
	GLuint simulationShader = 0;			//!< Compute shader

	GLuint bedrockBuffer = 0;				//!< Bedrock elevation buffer
	GLuint tempBedrockBuffer = 0;			//!< Output bedrock elevation buffer

	GLuint streamBuffer = 0;				//!< Water elevation buffer
	GLuint tempStreamBuffer = 0;				//!< Output water elevation buffer

	GLuint upliftBuffer = 0;		//!< Uplift buffer
	GLuint noiseBuffer = 0;			//!< Per-cell erodibility factor (k multiplier)

	int nx = 0;
	int ny = 0;
	int totalBufferSize = 0;
	int dispatchSize = 0;
	float ax = 0, ay = 0, bx = 0, by = 0;
	float cellDiagX = 0, cellDiagY = 0;
	std::vector<float> tmpData;
public:
	GPU_SPE() {};
	~GPU_SPE();

	void Init(const ScalarField2&);
	void ReloadShader();
	void Resize(const ScalarField2&);
	void Step(int n);
	void SetDt(float dt) const;

	void SetUplift(const ScalarField2& uplift) const;
	void SetNoise(const ScalarField2& noise) const;
	GLuint GetData() const;
	void GetData(ScalarField2& sf);
	void GetData(ScalarField2& sf, ScalarField2& sa);
};

#endif