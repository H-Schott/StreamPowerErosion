#pragma once

class Random
{
protected:
	int vec[48];
	int feed, tap, borrow;
	int V[571];
	int Y;

public:
	Random(int = 239);
	~Random();

	void Seed(int);
	int Integer(int);
	double Uniform();
	double Uniform(const double&);
	double Uniform(const double&, const double&);
public:
	static Random R239;
protected:
	int Int31();
	int SimpleInt31();
};
