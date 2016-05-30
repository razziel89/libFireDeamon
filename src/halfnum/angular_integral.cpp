/***********
This file is part of libFireDeamon.

Copyright (C) 2015,2016 by Torsten Sachse

libFireDeamon is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

libFireDeamon is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with libFireDeamon.  If not, see <http://www.gnu.org/licenses/>.
***********/
#include <cstdio>
#include <algorithm>
#include <cstring>
#include <string>
#include <cmath>

const static double Pi  = acos(-1.0L);
const static double S00 = 1.0L/(2.0L*sqrt(Pi));

//computes n!
inline long int f_factorial(int n){
	long int fact;
	fact=1;
	for (long int i=n;i>=2;--i){
		fact*=i;
	}
	return fact;
}

//computes n!!
inline long int f_double_factorial(int n){
	long int fact;
	fact=1;
	for (long int i=n;i>=2;i-=2){
		fact*=i;
	}
	return fact;
}


//computes n over k
inline int f_binomial (int n, int k){
	return int(f_factorial(n)/(f_factorial(n-k)*f_factorial(k)));
}

//return true if i is odd
inline bool f_odd(int i){
	if (i%2 == 1 && i >= 0){return true;}
	else{return false;}
}

//return true if i is even
inline bool f_even(int i){
	if (i%2 == 0 && i >= 0){return true;}
	else{return false;}
}

//calculation of inner conditional factor
inline double f_u_coeff_inner(int m, int lx){
	double result;
	if 	(m>0  && f_even(abs(abs(m)-lx))){
		result=1.0;
	}
	else{if	(m==0 && f_even(lx)){
		result=1.0/(sqrt(2));
	}
	else{if	(m<0  && f_odd(abs(abs(m)-lx))){
		result=1.0;
	}
	else{
		result=0;
	}}}
	return result;
}

//calculation of inner factor, i.e. the sum over the product of 2 binomials
//with a -1 or 1, depending on the case
inline int f_u_coeff_middle(int j, int m, int lx){
	int result=0;
	for(int k=0;k<=j;++k){
		if (not(j<0 or lx-2*k<0 or j-k<0 or abs(m)-lx+2*k<0)){
			result+=f_binomial(j,k)*f_binomial(abs(m),lx-2*k)*pow(-1,floor((abs(m)-lx+2*k)/2.0));
		}
		else{
			result+=0;
		}
	}
	return result;
}

//the product of 2 binomials with the quotient of 2 factorials
inline double f_u_coeff_factorial_fraction_and_binomials(int l, int i, int j, int m){
	double result;
	if (2*l-2*i < 0 || l-abs(m)-2*i < 0 || l<0 || i<0 || j<0 || l-i<0 || i-j<0){
		result=0;
	}
	else{
		result=f_binomial(l,i)*f_binomial(i,j)*pow(-1,i)*f_factorial(2*l-2*i)/(1.0*f_factorial(l-abs(m)-2*i));
	}
	return result;
}

//this function evaluates equation 15 of the given publication
inline double f_equation_15(int i, int j, int k){
	long double result;
	if ( f_even(i) && i>=0 && f_even(j) && j>=0 && f_even(k) && k>=0 ){
		result=1.0*f_double_factorial(i-1)/f_double_factorial(i+j+k+1)*f_double_factorial(j-1)*f_double_factorial(k-1);
	}
	else{
		result=0;
	}
	return double(4*Pi*result);
}

//the top-level function putting the others together
inline double f_u_coeff (int l, int m, int lx, int ly, int lz){
	double result=0;
	int j=lx+ly-abs(m);
	if ( j%2==0 && j>=0 ){
		j/=2;
		double inner_factor=f_u_coeff_inner(m,lx);
		if (not(inner_factor==0)){
			double prefactor=sqrt((2*l+1)/(2*Pi)*f_factorial(l-abs(m))/f_factorial(l+abs(m)))/(pow(2,l)*f_factorial(l));
			double outer_sum=0;
			for (int i=j;2*i<=(l-abs(m));++i){
				double factorial_fraction_and_binomials=f_u_coeff_factorial_fraction_and_binomials(l,i,j,m);
				if (not(factorial_fraction_and_binomials==0)){
					outer_sum+=factorial_fraction_and_binomials*f_u_coeff_middle(j,m,lx);
				}
			}
			result=prefactor*outer_sum*inner_factor;
		}
	}
	return result;
}

//the top-level funtion for the v-coefficients
inline double f_v_coeff (int l, int m, int lx, int ly, int lz){
	double result=0;
	if (l<=lx+ly+lz){
		for (int i=0;i<=l;++i){
			for (int j=0;j<=l-i;++j){
				result+=f_u_coeff(l,m,i,j,l-i-j)*f_equation_15(lx+i,ly+j,lz+l-i-j);
			}
		}
	}
	return result;
}

////populate lxyz with (lx,ly,lz)
////this construct is used to map (lx,ly,lz) with sum((lx,ly,lz))==l to one index
//inline void f_populate_lxyz (int lmax, int lxyz[]){
//	int count_lxyz=0;
//	for (int l=0; l<=lmax; ++l){
//		for (int lx=0;lx<=l;++lx){
//			for (int ly=0;ly<=l;++ly){
//				for (int lz=0;lz<=l;++lz){
//					if (lx+ly+lz==l){
//						lxyz[0+count_lxyz*3]=lx;
//						lxyz[1+count_lxyz*3]=ly;
//						lxyz[2+count_lxyz*3]=lz;
//						++count_lxyz;
//					}
//				}
//			}
//		}
//	}
//	return;
//}
//
////populate alphaxyz with (alphax,alphay,alphaz)
////this construct is used to map (\alpha_x,\alpha_y,\alpha_z) with \alpha_i<=a_i and sum((ax,ay,az))==a to one index
//inline void f_populate_alphaxyz (int a, int alphaxyz[]){
//	int count_alphaxyz = 0;
//	for (int alphax=0;alphax<=a;++alphax){
//		for (int alphay=0;alphay<=a;++alphay){
//			for (int alphaz=0;alphaz<=a;++alphaz){
//				if (alphax+alphay+alphaz <= a){
//					alphaxyz[0+count_alphaxyz*3]=alphax;
//					alphaxyz[1+count_alphaxyz*3]=alphay;
//					alphaxyz[2+count_alphaxyz*3]=alphaz;
//					++count_alphaxyz;
//				}
//			}
//		}
//	}
//	return;
//}
//
////This does the same as the above functions with the only difference
////that it counts for a given ax, ay and az
//int f_count_alphaxyz_orbital (int ax, int ay, int az){
//	int nr_alphaxyz = 0;
//	for (int alphax=0;alphax<=ax;++alphax){
//		for (int alphay=0;alphay<=ay;++alphay){
//			for (int alphaz=0;alphaz<=az;++alphaz){
//				++nr_alphaxyz;
//			}
//		}
//	}
//	return nr_alphaxyz;
//}
//
////This does the same as the above functions with the only difference
////that it populates for a given ax, ay and az.
////This way, the appropriate alpha indeces can be found that are required to
////acces the appropriate omega-coefficients for a certain orbital type (given
////by ax, ay and az).
//void f_populate_alphaxyz_orbital (int a, int ax, int ay, int az, int* alphaxyz, int* a_minus_alphaxyz, int* counts, long long int* binomials){
//	int count_alphaxyz_orbital = 0;
//	int count_alphaxyz = 0;
//	for (int alphax=0;alphax<=a;++alphax){
//		for (int alphay=0;alphay<=a;++alphay){
//			for (int alphaz=0;alphaz<=a;++alphaz){
//				if (alphax <= ax && alphay <= ay && alphaz <= az){
//					alphaxyz[0+count_alphaxyz_orbital*3]=alphax;
//					alphaxyz[1+count_alphaxyz_orbital*3]=alphay;
//					alphaxyz[2+count_alphaxyz_orbital*3]=alphaz;
//					a_minus_alphaxyz[0+count_alphaxyz_orbital*3]=ax-alphax;
//					a_minus_alphaxyz[1+count_alphaxyz_orbital*3]=ay-alphay;
//					a_minus_alphaxyz[2+count_alphaxyz_orbital*3]=az-alphaz;
//					binomials[count_alphaxyz_orbital] = f_binomial(ax,alphax) * f_binomial(ay,alphay) * f_binomial(az,alphaz);
//					counts[count_alphaxyz_orbital]=count_alphaxyz;//these are the indices that
//					//associate the alphaxyz with the elements of the Omega coefficients
//					++count_alphaxyz_orbital;
//				}
//				if (alphax+alphay+alphaz <= a){
//					++count_alphaxyz;
//				}
//			}
//		}
//	}
//	return;
//}

//Return the square of an integer. The built-in function pow can
//do this only for doubles, but I want it for integers too.
//The actual problem with pow is that it appears to return a double
//even if not necessary.
inline int square (int i){
	return i*i;
}

#define LMAXP1 4
class AngInt{
    private:
        double* m_integrals;
    public:
        AngInt(){
            size_t max_nr_elements = LMAXP1 * LMAXP1 * LMAXP1 * (3*LMAXP1+1);
            m_integrals = (double*) malloc(max_nr_elements * sizeof(double));
            memset(m_integrals, 0.0, max_nr_elements * sizeof(double));

	        double u = f_u_coeff (0, 0, 0, 0, 0) / S00;

	        for (int i=0; i<LMAXP1; ++i){
	        	for (int j=0; j<LMAXP1; ++j){
	        		for (int k=0; k<LMAXP1; ++k){
                        int index = k + j*(LMAXP1) + i*(LMAXP1*LMAXP1);
                        if ( i+j+k < LMAXP1){
	        			    for (int lambda=0; lambda<i+j+k+1; ++lambda){
	        			    	for (int mu_plus_lambda=0; mu_plus_lambda<2*lambda+1; ++mu_plus_lambda){
                                    m_integrals[index + mu_plus_lambda+square(lambda)] = u * f_v_coeff (lambda, mu_plus_lambda-lambda, i, j, k);
	                            } //mu_plus_lambda
	                        } //lambda
                        } //check
	                } //k
	            } //j
	        } //i
        }
        double GetInt(unsigned int lambda, int mu, unsigned int i, unsigned int j, unsigned int k){
            if (i+j+k>=LMAXP1 || mu*mu>lambda*lambda){
                throw;
            }
            if (lambda>i+j+k){
                return 0.0;
            }else{
                return m_integrals[k + j*(LMAXP1) + i*(LMAXP1*LMAXP1) + mu+lambda+square(lambda)];
            }
        }
        ~AngInt(){
            free(m_integrals);
        }
};
#undef LMAXP1
