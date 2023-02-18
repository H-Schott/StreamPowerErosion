#include "spe_shader.h"


GPU_SPE::~GPU_SPE() {
	glDeleteBuffers(1, &bedrockBuffer);
	glDeleteBuffers(1, &tempBedrockBuffer);

	glDeleteBuffers(1, &streamBuffer);
	glDeleteBuffers(1, &tempStreamBuffer);

	glDeleteBuffers(1, &upliftBuffer);

	release_program(simulationShader);
}

void GPU_SPE::Init(const ScalarField2& hf) {
	// Prepare data for first step
	nx = hf.GetSizeX();
	ny = hf.GetSizeY();
	totalBufferSize = hf.VertexSize();
	dispatchSize = (max(nx, ny) / 8) + 1;

	tmpData.resize(totalBufferSize);
	for (int i = 0; i < totalBufferSize; i++)
		tmpData[i] = hf.at(i);

	// Prepare shader & Init buffer - Just done once
	std::string fullPath = "../Shaders/spe_shader.glsl";

	simulationShader = read_program(fullPath.c_str());

	if (bedrockBuffer == 0) glGenBuffers(1, &bedrockBuffer);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, bedrockBuffer);
	glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * totalBufferSize, &tmpData.front(), GL_STREAM_READ);

	if (tempBedrockBuffer == 0) glGenBuffers(1, &tempBedrockBuffer);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, tempBedrockBuffer);
	glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * totalBufferSize, &tmpData.front(), GL_STREAM_READ);

	if (streamBuffer == 0) glGenBuffers(1, &streamBuffer);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, streamBuffer);
	glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * totalBufferSize, &tmpData.front(), GL_STREAM_READ);

	if (tempStreamBuffer == 0) glGenBuffers(1, &tempStreamBuffer);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, tempStreamBuffer);
	glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * totalBufferSize, &tmpData.front(), GL_STREAM_READ);

	if (upliftBuffer == 0) glGenBuffers(1, &upliftBuffer);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, upliftBuffer);
	glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * totalBufferSize, &tmpData.front(), GL_STREAM_READ);

	// Uniforms - just once
	glUseProgram(simulationShader);

	Box2 box = hf.Array2::GetBox();
	Vector2 cellDiag = hf.CellDiagonal();
	glUniform1i(glGetUniformLocation(simulationShader, "nx"), nx);
	glUniform1i(glGetUniformLocation(simulationShader, "ny"), ny);
	glUniform2f(glGetUniformLocation(simulationShader, "cellDiag"), float(cellDiag[0]), float(cellDiag[1]));
	glUniform2f(glGetUniformLocation(simulationShader, "a"), float(box[0][0]), float(box[0][1]));
	glUniform2f(glGetUniformLocation(simulationShader, "b"), float(box[1][0]), float(box[1][1]));
	
	glUseProgram(0);
}

void GPU_SPE::Step(int n) {

}

void GPU_SPE::SetDt(float dt) const {

}

void GPU_SPE::SetUplift(const ScalarField2& uplift) const {

}

GLuint GPU_SPE::GetData() const {
	return bedrockBuffer;
}

