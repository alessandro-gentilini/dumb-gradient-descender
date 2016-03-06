#include <iostream>
#include <vector>
#include <cmath>

double euclidean_distance( const std::vector<double>& a, const std::vector<double>& b )
{
	double d = 0;
	for( size_t i = 0; i<a.size(); i++) {
		d += std::pow(a[i]-b[i],2);
	}
	return sqrt(d);
}

std::vector<double> operator-(const std::vector<double>& a,const std::vector<double>& b)
{
	std::vector<double> result(a.size());
	for( size_t i = 0; i<a.size(); i++) {
		result[i] = a[i]-b[i];
	}	
	return result;
}

std::vector<double> operator*(const double& a,const std::vector<double>& b)
{
	std::vector<double> result(b.size());
	for( size_t i = 0; i<b.size(); i++) {
		result[i] = a*b[i];
	}	
	return result;
}

template < typename FUN >
std::vector<double> gradient(const FUN& fun, const std::vector<double>& x, double delta)
{
	std::vector<double> result(x.size());
	for ( size_t i = 0; i<x.size(); i++) {
		std::vector<double> dx(x);
		dx[i] += delta;
		result[i] = (fun(dx)-fun(x))/delta;
	}
	return result;
}

std::ostream& operator<<(std::ostream& o, const std::vector<double>& v)
{
	o << v[0] << " " << v[1] << " " << v[2];
	return o;
}

double f(const std::vector<double>& x)
{
	return x[0]*x[0]+4*x[1]*x[1];
}

template < typename FUN >
std::vector<double> dumb_gradient_descent(const std::vector<double>& initial_estimate,const FUN& fun,const double& lambda,const double& delta,const double& epsilon)
{
	std::vector<double> theta_new;
	std::vector<double> theta_old( initial_estimate );
	std::cout << "v=[" << theta_old << "\n";
	bool stop = false;
	do {
		theta_new = theta_old - lambda * gradient(f, theta_old, delta);
		if ( euclidean_distance(theta_new,theta_old) < epsilon ) {
			stop = true;
		}
		theta_old = theta_new;
		std::cout << theta_old << "\n";
	} while (!stop);
	return theta_old;
}

#include <cassert>
void minitest()
{
	std::vector< double > theta_new={1,2,3};
	std::vector< double > theta_old={1,2,3};
	double epsilon = .001;	
	assert(euclidean_distance(theta_old,theta_new)<epsilon);
	std::vector<double> zero={0,0,0};
	assert(theta_old-theta_new==zero);
	std::vector<double> doubles={2,4,6};
	assert(2*theta_new==doubles);

	double delta = .001;
	auto g(gradient(f,theta_new,delta));
}

int main(int, char**)
{
	std::vector< double > theta_new;
	std::vector< double > theta_old={51,.1,0};
	double epsilon = .001;
	double lambda = .1;
	double delta = lambda/10;

	minitest();
	
	dumb_gradient_descent(theta_old,f,lambda,delta,epsilon);

	std::cout << "];\n";
	std::cout << "x = linspace(-2,2);\n";
	std::cout << "y = linspace(-2,2);\n";
	std::cout << "[X,Y] = meshgrid(x,y);\n";
	std::cout << "Z = X.^2+4*Y.^2;\n";
	std::cout << "contour(X,Y,Z)\n";
	std::cout << "hold on \n";
	std::cout << "plot(v(:,1),v(:,2))\n";

	return 0;
}