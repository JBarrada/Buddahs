#include "BMP.h"

#include <string>
#include <cstring>
#include <vector>
#include <cstdlib>

#include <stdint.h>
#include <complex>
#include <cmath>

using namespace std;

const int bmp_width = 2000;
const int bmp_height = 2000;
const double bmp_width_2 = bmp_width / 2.0;
const double bmp_height_2 = bmp_height / 2.0;
const double aspect = bmp_width / (double)bmp_height;

const complex<double> center(0.0, 0.0);
const double zoom = 4.0;
const double exp_zoom = exp(zoom * log(1.1));

const complex<double> vp_center(-1.18, 0.0);
const double vp_zoom = -36.0;
const double vp_exp_zoom = exp(vp_zoom * log(1.1));

const int max_calc = 1000;

double image_buffer[bmp_width][bmp_height][3];
double image_max[3] = {0,0,0};

BMP mandelbrot(bmp_width, bmp_height);
BMP buddahbrot(bmp_width, bmp_height);


vector <complex<double> > z_samples;
vector <double> c_samples; // contributions

const int metro_threads = 200;
double l[metro_threads];
double o[metro_threads];

complex<double> current_orbit[max_calc];
int current_orbit_length = 0;

double complex_distance(complex<double> a, complex<double> b) {
	return ((a.real() - b.real())*(a.real() - b.real()) + (a.imag() - b.imag())*(a.imag() - b.imag()));
}

double random_range(double low, double high) {
	double s = fabs(low - high);
	double f = rand() / (double)RAND_MAX;
	return low + (s * f);
}

void random_complex_range(complex<double> &c, double r) {
	while (1) {
		c = *(new complex<double>(random_range(-r, r), random_range(-r, r)));
		
		if ((c.real()*c.real() + (c.imag()*c.imag())) < (r*r)) {
			return;	
		}
	}
}

void random_complex(complex<double> &c) {
	random_complex_range(c, 2);
}

void complex_to_screen(complex<double> value, int &x, int &y) {
	x = (((value.real() - center.real()) / exp_zoom) + 1.0) * bmp_width_2;
	y = ((((value.imag() * aspect) - center.imag()) / exp_zoom) + 1.0) * bmp_height_2;
}

void complex_to_screen_vp(complex<double> value, int &x, int &y) {
	x = (((value.real() - vp_center.real()) / vp_exp_zoom) + 1.0) * bmp_width_2;
	y = ((((value.imag() * aspect) - vp_center.imag()) / vp_exp_zoom) + 1.0) * bmp_height_2;
}

bool on_screen(int x, int y) {
	return ((x >= 0 && x < bmp_width) && (y >= 0 && y < bmp_height));
}

bool evaluate(complex<double> z, int max_calc_local) {
	current_orbit_length = 0;
	
	complex<double> zc(z);
	for (int i = 0; i < max_calc_local; i++) {
		z = (z * z) + zc;
		if ((z.real()*z.real() + z.imag()*z.imag()) > 4.0) {
			return true;
		}
		
		current_orbit[current_orbit_length++] = complex<double>(z);
		
		if (current_orbit_length >= max_calc) {
			return false;
		}
	}
	return false;
}

complex<double> mutate(complex<double> &c) {
	complex<double> n(c);
	
	if (random_range(0, 5) < 4) {
		double r1 = (1.0 / exp_zoom) * 0.0001;
		double r2 = (1.0 / exp_zoom) * 0.1;
		double phi = random_range(0, 1) * 2.0 * 3.1415926;
		double r = r2 * exp(-log(r2 / r1) * random_range(0, 1));
		
		n += complex<double>(r / cos(phi), r / sin(phi));
		return n;
	} else {
		random_complex(n);
		return n;
	}
}

double get_contribution() {
	int x, y;
	
	double contribution = 0;

	for (int i = 0; i < current_orbit_length; i++) {
		complex_to_screen_vp(current_orbit[i], x, y);
		
		if (on_screen(x, y)) {
			contribution++;
		}
	}

	if (current_orbit_length != 0) {
		return contribution / (double)(current_orbit_length);
	} else {
		return 0;
	}
}

double get_transition_probability(double q1, double q2, double olen1, double olen2) {
	return (1.0-(q1-olen1)/q1) / (1.0-(q2-olen2)/q2);
}

bool find_initial_sample(complex<double> &c, double x, double y, double rad, int f) {
	//printf("find_initial_sample(c=(%0.2f, %0.2f), x=%0.2f, y=%0.2f, rad=%0.2f, f=%d)\n", c.real(), c.imag(), x, y, rad, f);
	if (f > 500)
		return false;
	
	complex<double> tmp, seed;
	
	double closest = 1e20;
	
	for (int i = 0; i < 200; i++) {
		random_complex_range(tmp, rad);
		tmp += *(new complex<double>(x, y));
		
		if (!evaluate(tmp, max_calc)) {
			continue;
		}
		
		if (get_contribution() > 0.0) {
			c = tmp;
			return true;
		}
		
		for (int q = 0; q < current_orbit_length; q++) {
			double d = complex_distance(current_orbit[q], center);
			if (d < closest) {
				closest = d;
				seed = tmp;
			}
		}
	}
	return find_initial_sample(c, seed.real(), seed.imag(), rad / 2.0, f+1);
}

void build_initial_sample_points() {
	z_samples.clear();
	c_samples.clear();
	
	for (int i = 0; i < metro_threads; i++) {
		complex<double> m(0.0, 0.0);
		
		l[i] = o[i] = 1;
		
		if (!find_initial_sample(m, 0, 0, 2.0, 0)) {
			printf("couldn't find seed %d\n", i);
			continue;
		}

		evaluate(m, max_calc);
		z_samples.push_back(m);
		c_samples.push_back(get_contribution());
	}
}

void warmup() {
	complex<double> n_sample;
	double n_contrib;
	
	for (int s = 0; s < z_samples.size(); s++) {
		for (int e = 0; e < 100; e++) {
			int calc_length = max_calc / 100;
			
			n_sample = mutate(z_samples[s]);
			
			if (!evaluate(n_sample, max_calc)) {
				continue;
			}
			
			n_contrib = get_contribution();
			
			if (n_contrib == 0) {
				continue;
			}
			
			double t1 = get_transition_probability(calc_length, l[s], current_orbit_length, o[s]);
			double t2 = get_transition_probability(l[s], calc_length, o[s], current_orbit_length);
			
			double alpha = min((double)1.0, exp(log(n_contrib*t1)-log(c_samples[s]*t2)));
			double r = random_range(0, 1.0);
			
			if (alpha > r) {
				c_samples[s] = n_contrib;
				z_samples[s] = n_sample;
				
				l[s] = calc_length;
				o[s] = current_orbit_length;
			}
		}
	}
}

void draw_orbit(int color, double g, complex<double> *orbit, int len) {
	int x, y;
	
	for (int i = 0; i < len; i++) {
		complex_to_screen_vp(orbit[i], x, y);
		if (on_screen(x, y)) {
			image_buffer[x][y][color] += g;
			
			if (image_buffer[x][y][color] > image_max[color]) {
				image_max[color] = image_buffer[x][y][color];
			}
		}
	}
}

void render() {
	complex<double> n_sample;
	double n_contrib;
	
	/*
	uint8_t pixel[3];
	pixel[0] = 255;
	pixel[1] = 0;
	pixel[2] = 255;
	*/
	
	int max_episodes = 1000000;
	
	for (int episode = 0; episode < max_episodes; episode++) {
		for (int s = 0; s < z_samples.size(); s++) {
			int calc_length = max_calc / 100;
			
			n_sample = mutate(z_samples[s]);
			
			//if (!evaluate(n_sample, i*pow(10, calc_length))) {
			if (!evaluate(n_sample, max_calc)) {
				continue;
			}
			
			n_contrib = get_contribution();
			
			if (n_contrib == 0) {
				continue;
			}
			
			double t1 = get_transition_probability(calc_length, l[s], current_orbit_length, o[s]);
			double t2 = get_transition_probability(l[s], calc_length, o[s], current_orbit_length);
			
			double alpha = min((double)1.0, exp(log(n_contrib*t1)-log(c_samples[s]*t2)));
			double r = random_range(0, 1.0);
			
			
			//printf("alpha=%f  r=%f\n", alpha, r);
			
			if (alpha > r) {
				c_samples[s] = n_contrib;
				z_samples[s] = n_sample;
				
				l[s] = calc_length;
				o[s] = current_orbit_length;
				
				
				int color = 2;
				if (current_orbit_length < (max_calc / 50.0))
					color = 0;
				if (current_orbit_length < (max_calc / 10.0))
					color = 1;
				
				draw_orbit(color, 1, current_orbit, current_orbit_length);
				
				/*
				int x, y;
				complex_to_screen(n_sample, x, y);
				mandelbrot.set_pixel(x, y, pixel);
				*/
			} else {
				int color = 2;
				if (current_orbit_length < (max_calc / 50.0))
					color = 0;
				if (current_orbit_length < (max_calc / 10.0))
					color = 1;
				
				draw_orbit(color, 1, current_orbit, current_orbit_length);
			}
		}
		
		if (episode % (max_episodes / 100) == 0) {
			printf("%0.1f%%\n", ((episode * 100.0) / (double)max_episodes));
		}
	}
}

void draw_mandelbrot() {
	uint8_t pixel[3];
	
	double x_norm, y_norm;
	
	for (int x = 0; x < bmp_width; x++) {
		for (int y = 0; y < bmp_height; y++) {
			x_norm = (x / bmp_width_2) - 1.0;
			y_norm = (y / bmp_height_2) - 1.0;
			
			complex<double> z(x_norm*exp_zoom + center.real(), (y_norm*exp_zoom + center.imag()) / aspect);
			
			if (evaluate(z, max_calc / 5)) {
				pixel[0] = pixel[1] = pixel[2] = 255;
				mandelbrot.set_pixel(x, y, pixel);
			}
		}
	}
	
	/*
	int x, y;
	for (int s = 0; s < z_samples.size(); s++) {
		complex_to_screen(z_samples[s], x, y);
		
		//printf("s=(%0.4f, %0.4f) : %0.4f\n", z_samples[s].real(), z_samples[s].imag(), c_samples[s]);
		
		if (on_screen(x, y)) {
			pixel[0] = 255;
			pixel[1] = 0;
			pixel[2] = 255;
			mandelbrot.set_pixel(x, y, pixel);
		}
	}
	*/
}

int main() {
	build_initial_sample_points();
	//warmup();
	render();
	
	//draw_mandelbrot();
	
	uint8_t pixel[3];
	double global_image_max = max(max(image_max[0], image_max[1]), image_max[2]);
	for (int x = 0; x < bmp_width; x++) {
		for (int y = 0; y < bmp_height; y++) {
			pixel[0] = (image_buffer[x][y][0] / image_max[0]) * 255;
			pixel[1] = (image_buffer[x][y][1] / image_max[1]) * 255;
			pixel[2] = (image_buffer[x][y][2] / image_max[2]) * 255;
			buddahbrot.set_pixel(x, y, pixel);
		}
	}
	printf("max %f\n", image_max);
	
	buddahbrot.save("buddha.bmp");
	

	//draw_mandelbrot();
	
	mandelbrot.save("test.bmp");
	
	return 0;
}