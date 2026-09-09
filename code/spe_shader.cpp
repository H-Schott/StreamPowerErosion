#include "spe_shader.h"


GPU_SPE::~GPU_SPE() {
	glDeleteBuffers(1, &bedrockBuffer);
	glDeleteBuffers(1, &tempBedrockBuffer);

	glDeleteBuffers(1, &streamBuffer);
	glDeleteBuffers(1, &tempStreamBuffer);

	glDeleteBuffers(1, &upliftBuffer);
	glDeleteBuffers(1, &noiseBuffer);

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

	std::vector<float> tmpZeros(totalBufferSize, 0.);

	// Prepare shader & Init buffer - Just done once
	std::string fullPath = std::string(PATH_TO_SRC_DIRECTORY) + "data/shaders/spe_shader.glsl";

	simulationShader = read_program(fullPath.c_str());

	if (bedrockBuffer == 0) glGenBuffers(1, &bedrockBuffer);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, bedrockBuffer);
	glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * totalBufferSize, &tmpData.front(), GL_STREAM_READ);

	if (tempBedrockBuffer == 0) glGenBuffers(1, &tempBedrockBuffer);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, tempBedrockBuffer);
	glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * totalBufferSize, &tmpZeros.front(), GL_STREAM_READ);

	if (streamBuffer == 0) glGenBuffers(1, &streamBuffer);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, streamBuffer);
	glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * totalBufferSize, &tmpZeros.front(), GL_STREAM_READ);

	if (tempStreamBuffer == 0) glGenBuffers(1, &tempStreamBuffer);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, tempStreamBuffer);
	glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * totalBufferSize, &tmpZeros.front(), GL_STREAM_READ);

	if (upliftBuffer == 0) glGenBuffers(1, &upliftBuffer);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, upliftBuffer);
	glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * totalBufferSize, &tmpZeros.front(), GL_STREAM_READ);

	std::vector<float> tmpOnes(totalBufferSize, 1.f);
	if (noiseBuffer == 0) glGenBuffers(1, &noiseBuffer);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, noiseBuffer);
	glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * totalBufferSize, &tmpOnes.front(), GL_STREAM_READ);

	// Uniforms - just once
	glUseProgram(simulationShader);

	Box2 box = hf.Array2::GetBox();
	Vector2 cellDiag = hf.CellDiagonal();
	ax = float(box[0][0]); ay = float(box[0][1]);
	bx = float(box[1][0]); by = float(box[1][1]);
	cellDiagX = float(cellDiag[0]); cellDiagY = float(cellDiag[1]);
	std::cout << cellDiagX << " " << cellDiagY << std::endl;
	glUniform1i(glGetUniformLocation(simulationShader, "nx"), nx);
	glUniform1i(glGetUniformLocation(simulationShader, "ny"), ny);
	glUniform2f(glGetUniformLocation(simulationShader, "cellDiag"), cellDiagX, cellDiagY);
	glUniform2f(glGetUniformLocation(simulationShader, "a"), ax, ay);
	glUniform2f(glGetUniformLocation(simulationShader, "b"), bx, by);
	
	glUseProgram(0);
}

void GPU_SPE::ReloadShader() {
	release_program(simulationShader);
	simulationShader = read_program((std::string(PATH_TO_SRC_DIRECTORY) + "data/shaders/spe_shader.glsl").c_str());

	glUseProgram(simulationShader);
	glUniform1i(glGetUniformLocation(simulationShader, "nx"), nx);
	glUniform1i(glGetUniformLocation(simulationShader, "ny"), ny);
	glUniform2f(glGetUniformLocation(simulationShader, "cellDiag"), cellDiagX, cellDiagY);
	glUniform2f(glGetUniformLocation(simulationShader, "a"), ax, ay);
	glUniform2f(glGetUniformLocation(simulationShader, "b"), bx, by);
	glUseProgram(0);
}

void GPU_SPE::Resize(const ScalarField2& hf) {
	nx = hf.GetSizeX();
	ny = hf.GetSizeY();
	totalBufferSize = hf.VertexSize();
	dispatchSize = (max(nx, ny) / 8) + 1;

	tmpData.resize(totalBufferSize);
	for (int i = 0; i < totalBufferSize; i++)
		tmpData[i] = hf.at(i);

	std::vector<float> tmpZeros(totalBufferSize, 0.f);
	std::vector<float> tmpOnes(totalBufferSize, 1.f);

	glBindBuffer(GL_SHADER_STORAGE_BUFFER, bedrockBuffer);
	glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * totalBufferSize, tmpData.data(), GL_STREAM_READ);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, tempBedrockBuffer);
	glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * totalBufferSize, tmpZeros.data(), GL_STREAM_READ);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, streamBuffer);
	glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * totalBufferSize, tmpZeros.data(), GL_STREAM_READ);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, tempStreamBuffer);
	glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * totalBufferSize, tmpZeros.data(), GL_STREAM_READ);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, upliftBuffer);
	glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * totalBufferSize, tmpZeros.data(), GL_STREAM_READ);
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, noiseBuffer);
	glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * totalBufferSize, tmpOnes.data(), GL_STREAM_READ);

	Vector2 cellDiag = hf.CellDiagonal();
	cellDiagX = float(cellDiag[0]); cellDiagY = float(cellDiag[1]);

	glUseProgram(simulationShader);
	glUniform1i(glGetUniformLocation(simulationShader, "nx"), nx);
	glUniform1i(glGetUniformLocation(simulationShader, "ny"), ny);
	glUniform2f(glGetUniformLocation(simulationShader, "cellDiag"), cellDiagX, cellDiagY);
	glUseProgram(0);
}

void GPU_SPE::Step(int n) {

	for (int i = 0; i < n; i++) {

		glUseProgram(simulationShader);
		glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, bedrockBuffer);
		glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, streamBuffer);
		glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, tempBedrockBuffer);
		glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, tempStreamBuffer);
		glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 4, upliftBuffer);
		glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 5, noiseBuffer);

		glDispatchCompute(dispatchSize, dispatchSize, 1);
		glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

		// dual buffering
		std::swap(bedrockBuffer, tempBedrockBuffer);
		std::swap(streamBuffer, tempStreamBuffer);
	}

	glUseProgram(0);
}

void GPU_SPE::SetDt(float dt) const {
	glUseProgram(simulationShader);
	glUniform1f(glGetUniformLocation(simulationShader, "dt"), dt);
	glUseProgram(0);
}

void GPU_SPE::SetUplift(const ScalarField2& uplift) const {
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, upliftBuffer);
	glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * uplift.VertexSize(), &uplift.GetFloatData()[0], GL_STREAM_READ);
}

void GPU_SPE::SetNoise(const ScalarField2& noise) const {
	glBindBuffer(GL_SHADER_STORAGE_BUFFER, noiseBuffer);
	glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * noise.VertexSize(), &noise.GetFloatData()[0], GL_STREAM_READ);
}

GLuint GPU_SPE::GetData() const {
	return bedrockBuffer;
}

void GPU_SPE::GetData(ScalarField2& sf) {
	glGetNamedBufferSubData(bedrockBuffer, 0, sizeof(float) * totalBufferSize, tmpData.data());

	for (int i = 0; i < totalBufferSize; i++)
		sf[i] = double(tmpData[i]);

	/*double low, high;
	sf.GetRange(low, high);
	std::cout << low << " " << high << std::endl;*/
}

void GPU_SPE::GetData(ScalarField2& sf, ScalarField2& sa) {
	glGetNamedBufferSubData(bedrockBuffer, 0, sizeof(float) * totalBufferSize, tmpData.data());

	for (int i = 0; i < totalBufferSize; i++)
		sf[i] = double(tmpData[i]);

	glGetNamedBufferSubData(streamBuffer, 0, sizeof(float) * totalBufferSize, tmpData.data());

	for (int i = 0; i < totalBufferSize; i++)
		sa[i] = double(tmpData[i]);

	/*double low, high;
	sa.GetRange(low, high);
	std::cout << low << " " << high << std::endl;*/
}


