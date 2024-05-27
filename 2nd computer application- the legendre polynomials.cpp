#include <stdio.h>
#include <math.h>
#define M_PI 3.14159265358979323846

double P1(int n, double x){ //define a function for Pn(x)to calculate the Legendre polynomial using the recursive formula of the Legendre polynomial.takes two parameters n and x
	if (n == 0){
		return 1;
	}
	if (n == 1){
		return x;
	}
	if (n >= 2){
    	double p0 = 1; // Initialize p0 and p1 to the first two Legendre polynomials
        double p1 = x;
        double p;
        for (int i = 2; i <= n; i++) { //for loop to iterate over the values of i from 2 to n and calculates the polynomial using the formula of the Legendre polynomial.
            p = ((2.0 * i - 1.0) / i) * x * p1 - ((i - 1.0) / i) * p0;
            p0 = p1;
            p1 = p;
        }
        return p;
	}
}

double P2(int n, double x){  // define a function for Pn(x) to calculate the Legendre polynomial using the asymptotic expansion,takes two parameters n and x
    double arg = (n + 0.5) * acos(x) - M_PI / 4; // Compute the argument of the cosine function
    if (fabs(x) <= 1.0) {
	    return sqrt(2/(M_PI*n*sin(acos(x))))*cos(arg); // Compute the Legendre polynomial using the asymptotic expansion
    }
}

 
typedef struct {
    double x;
    double y;
    double z;
} Vector;

double term(Vector r, Vector r_prime, double Q, int n, double x) { //define a function to compute the a-th term of the sum in the force calculation
    int a = 0;
    x = 0.8; // set x to 0.8
    double sum = 0;
    double magnitude1=sqrt(pow(r_prime.x,2)+pow(r_prime.y,2)+pow(r_prime.z,2)); // Magnitude of r_prime
    double magnitude2=sqrt(pow(r.x, 2) + pow(r.y, 2) + pow(r.z, 2)); // Magnitude of r 
    double dot_product = r.x * r_prime.x + r.y * r_prime.y + r.z * r_prime.z;
    double cos_theta= dot_product / (magnitude1 * magnitude2);  // Cosine of the angle between r and r_prime
    double theta = acos(cos_theta); // Angle between r and r_prime
    while (a != n){ // Compute the a-th term of the sum
        sum += (pow(magnitude1,a) / pow(magnitude2,(a + 1))) * P1(a,cos_theta);
        a++;
    }
    return sum;
} 


void computeForce(Vector r, Vector r_prime, double Q) { // define a function to compute the force between two point charges Q at r and Q' at r'
    double inverse_difference; 
    double magnitude1=sqrt(pow(r_prime.x,2)+pow(r_prime.y,2)+pow(r_prime.z,2)); // Magnitude of r_prime
    double magnitude2=sqrt(pow(r.x, 2) + pow(r.y, 2) + pow(r.z, 2)); // Magnitude of r

    // calculate the magnitude of the displacement vector r-r'
    double magnitude = sqrt(pow(r.x - r_prime.x, 2) + pow(r.y - r_prime.y, 2) + pow(r.z - r_prime.z, 2)); // calculate the magnitude of the displacement vector r-r'
    if (magnitude2> magnitude1){  //check for |r|>|r'|
        inverse_difference = 1 / (magnitude); //calculate the inverse distance of two vectors (1/(|r-r'|))
    }
    double x = 0.8; //set x to 0.8
    double calc_1 = inverse_difference;
    double dot_product = r.x * r_prime.x + r.y * r_prime.y + r.z * r_prime.z; // calculate the dot product of vectors to find the cos(Q)
    double cos_Q = dot_product / (magnitude1 * magnitude2); //evaluate cos(Q) 
    Q = acos(cos_Q);

    printf("The angle between the vectors is : %lf radian ", Q);
    printf("Cosine value of Q is : %lf\n", cos_Q);
    printf("1/|r-r'| = %lf\n", inverse_difference);
    double term_0 = (pow(magnitude1,0) / pow(magnitude2,(1))) * P1(0,cos_Q);
    //printf("The result for n = %d is %lf\n", 0, term_0);
    //printf("The result for n = %d is %lf\n", 1, general_term);
    double term_1 = (pow(magnitude1,1) / pow(magnitude2,(2))) * P1(1,cos_Q);
    int n = 2;
    double general_term = term(r, r_prime, Q, n, x);
    double error_2 = fabs(inverse_difference - general_term); //find the absolute error between the calculations 1/(|r-r'|) and |r|^n / |r|^(n+1)*Pn(cosQ)
    while(error_2 > 0.0001 ){
        general_term = term(r,r_prime, Q, n+1, x);
        error_2 = fabs(inverse_difference- general_term);
        printf("The result for n = %d is %lf\n", n, general_term);
        n++;
    }
    printf("The appropritae value of n is: %d\n", n-1);  //print the appropraite value for n that makes the absolute error less than 0.0001
}


int main() { 

double x;
printf("Enter an x value between -1 and 1: "); //taking input from user for x 
scanf("%lf",&x);

double pol[10]; //define an array to calculate the sum of the first 10 terms by using Pn(x) the Legendre polynomial calculation 
double sum = 0;
for (int i = 0; i < 10; i++){
	pol[i] = P1(i,x);
	//printf("%lf\n",pol[i]);
	sum += pol[i];
}

printf("The sum of first 10 terms=%lf\n", sum);
double res_1 = P1(50,x);
double res_2 = P2(50,x);
printf("Result of legendre polynomial calculation for n=50: %lf\n",P1(50,x));
printf("Result of the asymptotic expansion for n=50: %lf\n",P2(50,x));

int n = 1;

double error = fabs(res_1 - res_2); //calculate the error for the while loop till the error became less than 0.0001
while (error >= 0.0001) {
    n++;
    res_1 = P1(n, x);
    res_2 = P2(n, x);
    error = fabs(res_1 - res_2);
}
printf("The n value that makes the error less then 10^-4 = %d\n", n);
double wanted_n1 = P1(n,x);
double wanted_n2= P2(n,x);
printf("Result of legendre polynomial calculation for n=51: %lf\n",P1(51,x));
printf("Result of the asymptotic expansion for n=51: %lf\n",P2(51,x));

printf("Absolute error for n=51:%lf\n", error);

double Q;
double a1 = 3.2; //define the values for the vectors from the given file
double b1 = 4.0;
double c1 = 7.5;
double a2 = -2.1;
double b2 = 3.3;
double c2 = -4.9;
int a = 1;
// define r vector
Vector r = {a1, b1, c1}; // define r vector

// define r' vector
Vector r_prime = {a2, b2, c2}; // define r' vector

// print the coordinates of the vectors
//printf("r = (%f, %f, %f)\n", r.x, r.y, r.z);
//printf("r' = (%f, %f, %f)\n", r_prime.x, r_prime.y, r_prime.z);

computeForce(r, r_prime, Q); //call the function computeForce

return 0; //program exits code properly 
}

