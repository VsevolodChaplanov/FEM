#ifndef __CUBE_ELEMENT_CPP__
#define __CUBE_ELEMENT_CPP__

#include "CubeElem.h"

double CubeElem::phi1(const double &_xi, const double &_eta, const double &_zeta)
{
	return (1-_xi) * (1-_zeta) * (1-_eta);
}

double CubeElem::phi2(const double &_xi, const double &_eta, const double &_zeta)
{
	return _xi * (1-_eta) * (1-_zeta);
}

double CubeElem::phi3(const double &_xi, const double &_eta, const double &_zeta)
{
	return _xi * _eta * (1-_zeta);
}

double CubeElem::phi4(const double &_xi, const double &_eta, const double &_zeta)
{
	return (1-_xi) * _eta * (1-_zeta);
}

double CubeElem::phi5(const double &_xi, const double &_eta, const double &_zeta)
{
	return (1-_xi) * (1-_eta) * _zeta ;
}

double CubeElem::phi6(const double &_xi, const double &_eta, const double &_zeta)
{
	return _xi * (1-_eta) * _zeta ;
}

double CubeElem::phi7(const double &_xi, const double &_eta, const double &_zeta)
{
	return _xi * _eta * _zeta;
}

double CubeElem::phi8(const double &_xi, const double &_eta, const double &_zeta)
{
	return (1-_xi) * _eta * _zeta;
}

#endif