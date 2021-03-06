#include <iostream>
#include <vector>
#include <cmath>

double euclidean_distance( const std::vector<double>& a, const std::vector<double>& b )
{
    double d = 0;
    for ( size_t i = 0; i < a.size(); i++) {
        d += std::pow(a[i] - b[i], 2);
    }
    return sqrt(d);
}

double magnitude( const std::vector<double>& a )
{
    return euclidean_distance(a, {0., 0., 0.});
}

std::vector<double> operator-(const std::vector<double>& a, const std::vector<double>& b)
{
    std::vector<double> result(a.size());
    for ( size_t i = 0; i < a.size(); i++) {
        result[i] = a[i] - b[i];
    }
    return result;
}

std::vector<double> operator*(const double& a, const std::vector<double>& b)
{
    std::vector<double> result(b.size());
    for ( size_t i = 0; i < b.size(); i++) {
        result[i] = a * b[i];
    }
    return result;
}

template < typename FUN >
std::vector<double> gradient(const FUN& fun, const std::vector<double>& x, double delta)
{
    std::vector<double> result(x.size());
    for ( size_t i = 0; i < x.size(); i++) {
        std::vector<double> dx(x);
        dx[i] += delta;
        result[i] = (fun(dx) - fun(x)) / delta;
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
    return std::pow(x[0] - 1, 2) + 4 * std::pow(x[1] + 2, 2);
}

template < typename FUN >
std::vector<double> dumb_gradient_descent(const std::vector<double>& initial_estimate, const FUN& fun, const double& initial_lambda, const double& delta, const double& epsilon, const char* name)
{
    double lambda(initial_lambda);
    std::vector<double> theta_new;
    std::vector<double> theta_old( initial_estimate );
    std::cout << name << "=[" << theta_old << " " << fun(theta_old) << "\n";
    bool stop = false;
    bool decreasing;
    double g;
    do {
        auto g = gradient(fun, theta_old, delta);
        theta_new = theta_old - (lambda / magnitude(g)) * g;
        theta_new[1] = .5;
        if ( fun(theta_new) < fun(theta_old) ) {
            lambda *= 2;
            decreasing = true;
        } else {
            lambda *= .5;
            decreasing = false;
        }
        if ( euclidean_distance(theta_new, theta_old) < epsilon ) {
        // if ( magnitude(g) < epsilon ) {
            stop = true;
        }
        if ( decreasing ) {
            theta_old = theta_new;
        }
        std::cout << theta_old << " " << fun(theta_old) << "\n";
    } while (!stop);
    std::cout << "];";
    return theta_old;
}

#include <cassert>
void minitest()
{
    std::vector< double > theta_new = {1, 2, 3};
    std::vector< double > theta_old = {1, 2, 3};
    double epsilon = .001;
    assert(euclidean_distance(theta_old, theta_new) < epsilon);
    std::vector<double> zero = {0, 0, 0};
    assert(theta_old - theta_new == zero);
    std::vector<double> doubles = {2, 4, 6};
    assert(2 * theta_new == doubles);

    double delta = .001;
    auto g(gradient(f, theta_new, delta));
}

// http://git.savannah.gnu.org/cgit/gsl.git/tree/poly/solve_cubic.c
/* solve_cubic.c - finds the real roots of x^3 + a x^2 + b x + c = 0 */
#define SWAP(a,b) do { double tmp = b ; b = a ; a = tmp ; } while(0)

int
gsl_poly_solve_cubic (double a, double b, double c,
                      double *x0, double *x1, double *x2)
{
    double q = (a * a - 3 * b);
    double r = (2 * a * a * a - 9 * a * b + 27 * c);

    double Q = q / 9;
    double R = r / 54;

    double Q3 = Q * Q * Q;
    double R2 = R * R;

    double CR2 = 729 * r * r;
    double CQ3 = 2916 * q * q * q;

    if (R == 0 && Q == 0)
    {
        *x0 = - a / 3 ;
        *x1 = - a / 3 ;
        *x2 = - a / 3 ;
        return 3 ;
    }
    else if (CR2 == CQ3)
    {
        /* this test is actually R2 == Q3, written in a form suitable
           for exact computation with integers */

        /* Due to finite precision some double roots may be missed, and
           considered to be a pair of complex roots z = x +/- epsilon i
           close to the real axis. */

        double sqrtQ = sqrt (Q);

        if (R > 0)
        {
            *x0 = -2 * sqrtQ  - a / 3;
            *x1 = sqrtQ - a / 3;
            *x2 = sqrtQ - a / 3;
        }
        else
        {
            *x0 = - sqrtQ  - a / 3;
            *x1 = - sqrtQ - a / 3;
            *x2 = 2 * sqrtQ - a / 3;
        }
        return 3 ;
    }
    else if (R2 < Q3)
    {
        double sgnR = (R >= 0 ? 1 : -1);
        double ratio = sgnR * sqrt (R2 / Q3);
        double theta = acos (ratio);
        double norm = -2 * sqrt (Q);
        *x0 = norm * cos (theta / 3) - a / 3;
        *x1 = norm * cos ((theta + 2.0 * M_PI) / 3) - a / 3;
        *x2 = norm * cos ((theta - 2.0 * M_PI) / 3) - a / 3;

        /* Sort *x0, *x1, *x2 into increasing order */

        if (*x0 > *x1)
            SWAP(*x0, *x1) ;

        if (*x1 > *x2)
        {
            SWAP(*x1, *x2) ;

            if (*x0 > *x1)
                SWAP(*x0, *x1) ;
        }

        return 3;
    }
    else
    {
        double sgnR = (R >= 0 ? 1 : -1);
        double A = -sgnR * pow (fabs (R) + sqrt (R2 - Q3), 1.0 / 3.0);
        double B = Q / A ;
        *x0 = A + B - a / 3;
        return 1;
    }
}

double parabola_point_distance(double a, double b, double c, double x, double y)
{
    double t3 = 4 * std::pow(a, 2);
    double t2 = 6 * a * b;
    double t1 = 4 * a * c - 4 * a * y + 2 * std::pow(b, 2) + 2;
    double t0 = 2 * b * c - 2 * b * y - 2 * x;
    double x0, x1, x2;
    int n_roots = gsl_poly_solve_cubic(t2 / t3, t1 / t3, t0 / t3, &x0, &x1, &x2);
    double d = 0;
    switch (n_roots) {
    case 1:
        d = euclidean_distance({x, y}, {x0, a * std::pow(x0, 2) + b * x0 + c});
        break;
    case 3:
        double d1 = euclidean_distance({x, y}, {x0, a * std::pow(x0, 2) + b * x0 + c});
        double d2 = euclidean_distance({x, y}, {x1, a * std::pow(x1, 2) + b * x1 + c});
        double d3 = euclidean_distance({x, y}, {x2, a * std::pow(x2, 2) + b * x2 + c});
        d = std::min(d1, d2);
        d = std::min(d, d3);
        break;
    }
    return d;
}

class Cost
{
public:
    Cost( const std::vector<double>& x, const std::vector<double>& y )
        : x_(x), y_(y) {}

    double operator()(const std::vector<double>& theta) const
    {
        double acc = 0;
        for ( size_t i = 0; i < x_.size(); i++ ) {
            acc += parabola_point_distance(theta[0], theta[1], theta[2], x_[i], y_[i]);
        }
        return acc;
    }

    std::vector<double> x_;
    std::vector<double> y_;
};

int main(int, char**)
{
    std::vector< double > theta_old = { 1, 2, 3};
    double epsilon = .00001;
    double lambda = 10;
    double delta = .0001;

    minitest();

    //dumb_gradient_descent(theta_old,f,lambda,delta,epsilon);

    double a = .5;
    double b = .5;
    double c = .5;

    std::vector<double> x, y;

    for (int i = -10; i < 10; i++ ) {
        x.push_back(i/*+static_cast<double>(rand())/RAND_MAX*/);
        y.push_back(a * i * i + b * i + c/*+static_cast<double>(rand())/RAND_MAX*/);
    }


    Cost cost(x, y);

    std::vector< double > aa;
    std::vector< double > bb;
    std::vector< double > cc;

    size_t N = 1000;
    double step = 6. / N;
    for ( size_t i = 0; i < N; i++) {
        aa.push_back(-3 + i * step);
    }

    bb = aa;
    cc = aa;


    std::cout << "s=[";
    for ( size_t r = 0; r < aa.size(); r++ ) {
//        for ( size_t c = 0; c < bb.size(); c++ ) {
        for ( size_t p = 0; p < cc.size(); p++ ) {
            std::cout << aa[r] << " " << cc[p] << " " << cost({aa[r], .5, cc[p]}) << "\n";
        }
        //}
        std::cout << "\n";
    }
    std::cout << "];\n";

    std::string n1("v1");
    dumb_gradient_descent(theta_old, cost, lambda, delta, epsilon, n1.c_str());

    theta_old = { 1, 2, 1};
    std::string n2("v2");
    dumb_gradient_descent(theta_old, cost, lambda, delta, epsilon, n2.c_str());

    std::cout << "Z = reshape(s(:,3),[]," << N << ");\n";

    std::cout << "X = reshape(s(:,1),[]," << N << ");\n";

    std::cout << "Y = reshape(s(:,2),[]," << N << ");\n";


    std::cout << "contour(X,Y,Z,100)\n";
    std::cout << "xlabel('a')\n";
    std::cout << "ylabel('c')\n";
    std::cout << "hold on\n";
    std::cout << "plot(" << n1 << "(:,1)," << n1 << "(:,3),'k')\n";
    std::cout << "plot(" << n2 << "(:,1)," << n2 << "(:,3),'k')\n";
    std::cout << "plot(" << a << "," << c << ",'k*')\n";
    std::cout << "axis equal    \n";

    return 0;
}
