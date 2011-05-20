#ifndef __Rotation3D_h_
#define __Rotation3D_h_

class Rotation3D
{
public:
	Rotation3D();
	~Rotation3D();
	
	void UpdateMatrix();
	void UpdateParameter();
	double Get(int i);
	void Put(int i, double x);
	void Rotate(double& x, double& y, double& z);
	void SetCenter(double x,double y, double z);

private:
	double _rx;
	double _ry;
	double _rz;
	double _tx;
	double _ty;
	double _tz;
	double _matrix[4][4];
	double m_centerx;
	double m_centery;
	double m_centerz;
};


#endif // __Rotation3D_h_
