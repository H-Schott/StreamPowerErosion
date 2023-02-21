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

	int nx = 0;
	int ny = 0;
	int totalBufferSize = 0;					//!< Total buffer size defined as nx * ny
	int dispatchSize = 0;						//!< Single dispatch size
	std::vector<float> tmpData;			//!< Temporary array for retreiving GPU data
public:
	GPU_SPE() {};
	~GPU_SPE();

	void Init(const ScalarField2&);
	void Step(int n);
	void SetDt(float dt) const;

	void SetUplift(const ScalarField2& uplift) const;
	GLuint GetData() const;
	void GetData(ScalarField2& sf);
	void GetData(ScalarField2& sf, ScalarField2& sa);
};

#endif