/*

 g++ -Ofast -march=native -Wall -Wextra -Wpedantic -std=c++11 -o main main.cpp -lnetcdf_c++ /home/n/Dropbox/netcdfcxx4/netcdf-cxx4-4.2.1/cxx4/.libs/libnetcdf_c++4.so.1 -lnetcdf && ./main

 */

#include <iostream>
#include <functional>
#include <algorithm>
#include <tr1/functional>
#include <netcdf>
#include <string>
#include <bitset>
#include <tuple>
#include <cmath>

using namespace netCDF;
using namespace netCDF::exceptions;

static const int NLAT = 170;
static const int NLON = 180;
static const int NTIM = 24;

float data[NTIM][NLAT][NLON];

/* load data from tos_O1_2001-2002.nc to data */
bool load_nc();

static const int shl = (int) sizeof(size_t); /* simhash length */
static const int bs = 5; /* block size */
static const int binc = 1; /* block position increment */
static const int similarity_treshold = 10; /* how many bits can differ in simhash */
static const int search_window = 100;

/* stolen from boost */
void hash_combine(std::size_t& seed, const float& v) {
	std::tr1::hash<float> hasher;
	seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

void cyclic_shift(size_t &s) {
	size_t first_bit = s << 63;
	s >>= 1;
	s |= first_bit;
}

/******************************************************************************/

class Block {
public:
	int tim;
	int lat;
	int lon;
	std::size_t simhash;
	std::vector<Block> similar;

	Block(int t, int la, int lo, std::size_t s) :
			tim(t), lat(la), lon(lo), simhash(s) {
	}

	bool operator <(const Block& other) const {
		return (this->simhash < other.simhash);
	}

	bool operator==(const Block &other) const {
		/* this is probably not right. just compare pointers to instances */
		return (this->tim == other.tim && this->lat == other.lat
				&& this->lon == other.lon);
	}
};

/* TODO: something like BlockSimilarity class */

void block_diff(Block &b1, Block &b2) {
	int same = 0;
	int zero_same = 0;
	float diff = 0;

	for (int iti = 0; iti < bs; iti++) {
		for (int ilat = 0; ilat < bs; ilat++) {
			for (int ilon = 0; ilon < bs; ilon++) {
				float fb1 = data[b1.tim + iti][b1.lat + ilat][b1.lon + ilon];
				float fb2 = data[b2.tim + iti][b2.lat + ilat][b2.lon + ilon];

				if (fb1 == fb2) {
					if (fb1 == 0) {
						zero_same++;
						continue;
					}
					same++;
					continue;
				}

				diff += fabs(fb1 - fb2);
			}
		}
	}

	if (same > 0) {
		std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!! ";
	}

	std::bitset<64> simxor(b1.simhash ^ b2.simhash);
	std::cout << "zero_same=" << zero_same << ", same=" << same << ", diff="
			<< diff << ", xor=" << simxor.count() << std::endl;
}

void print_simhash(size_t simhash) {
	std::bitset<64> bs(simhash);
	std::cout << bs << std::endl;
}

void run_over_blocks() {

	std::tr1::hash<float> hasher;
	std::vector<Block> simhashes;
	std::vector<std::pair<Block, Block>> similar;

	for (int tim = 0; tim < NTIM - bs; tim += binc) {
		for (int lat = 0; lat < NLAT - bs; lat += binc) {
			for (int lon = 0; lon < NLON - bs; lon += binc) {
				std::vector<int> simhash_vec(64);
				int nonzero = 0;

				for (int iti = 0; iti < bs; iti++) {
					for (int ilat = 0; ilat < bs; ilat++) {
						for (int ilon = 0; ilon < bs; ilon++) {

							std::size_t seed = 0;
							float f = data[tim + iti][lat + ilat][lon + ilon];

							if (f == 0) {
								continue;
							}
							nonzero++;

							// https://stackoverflow.com/questions/19966041/getting-too-many-collisions-with-hash-combine
							hash_combine(seed, hasher(iti) * 2654435761);
							hash_combine(seed,
									hasher(ilat * 1000) * 2654435761);
							hash_combine(seed,
									hasher(ilon * 1000000) * 2654435761);
//							hash_combine(seed,
//									hasher(hasher(tim + lat * 1000 + lon * 100000))
//											* 2654435761);
							hash_combine(seed, hasher(f) * 2654435761);

							std::bitset<64> seed_bitset(seed);
//					std::cout << "seed_bitset: " << seed_bitset << std::endl;
							for (int j = 0; j < 64; j++) {
								if (seed_bitset[j]) {
									simhash_vec[j]++;
								} else {
									simhash_vec[j]--;
								}
							}

//					std::cout << "simhash_vec: ";
//					for (auto i = simhash_vec.begin(); i != simhash_vec.end(); ++i)
//					    std::cout << *i << ", ";
//					std::cout << std::endl;

						}
					}
				} /* for i 0..bs */

				/* if there was no nonzero value, don't save its simhash */
				if (nonzero == 0) {
					continue;
				}

				std::bitset<shl> simhash_bitset(0);
				for (int i = 0; i < 64; i++) {
					if (simhash_vec[i] >= 0) {
						simhash_bitset[i] = 1;
					} else {
						simhash_bitset[i] = 0;
					}
				}

				std::size_t simhash = simhash_bitset.to_ulong();
				simhashes.push_back(Block(tim, lat, lon, simhash));
			}
		}
	}

	std::cout << "now sort" << std::endl;

	for (int s = 0; s < 64; s++) {
		std::sort(simhashes.begin(), simhashes.end());

		for (int i = 0; i < (int) simhashes.size() - search_window; i++) {
			for (int sw = 0; sw < search_window; sw++) {
				std::bitset<64> simxor(
						simhashes[i].simhash ^ simhashes[i + sw].simhash);

				if (simxor.count() < similarity_treshold) {

					/* skip blocks too close together */
					/* TODO: make 5 a constant */
					if ((abs(simhashes[i].tim - simhashes[i + sw].tim) < 5)
							|| (abs(simhashes[i].lat - simhashes[i + sw].lat)
									< 5)
							|| (abs(simhashes[i].lon - simhashes[i + sw].lon)
									< 5)) {
						continue;
					}

					/* are they already similar? */
					if (std::find(simhashes[i].similar.begin(),
							simhashes[i].similar.end(), simhashes[i + sw])
							!= simhashes[i].similar.end()) {
						continue;
					}

					std::cout << "found almost same simhash" << std::endl;

					std::cout << "  first: ";
					std::cout << "tim=" << simhashes[i].tim << ", ";
					std::cout << "lon=" << simhashes[i].lon << ", ";
					std::cout << "lat=" << simhashes[i].lat;
					std::cout << std::endl;
					std::cout << "  ";
					print_simhash(simhashes[i].simhash);

					std::cout << "  second: ";
					std::cout << "tim=" << simhashes[i + sw].tim << ", ";
					std::cout << "lon=" << simhashes[i + sw].lon << ", ";
					std::cout << "lat=" << simhashes[i + sw].lat;
					std::cout << std::endl;
					std::cout << "  ";
					print_simhash(simhashes[i + sw].simhash);

					std::cout << "  ";
					block_diff(simhashes[i], simhashes[i + sw]);
					std::cout << std::endl;

					simhashes[i].similar.push_back(simhashes[i + sw]);
					simhashes[i + sw].similar.push_back(simhashes[i]);

					similar.push_back(
							std::make_pair(simhashes[i], simhashes[i + sw]));
				}
			}
		}

		for (int i = 0; i < (int) simhashes.size(); i++) {
			cyclic_shift(simhashes[i].simhash);
		}
	}

	std::cout << "found similar: " << similar.size() << std::endl;
}

/******************************************************************************/

int main() {
	if (!load_nc()) {
		return (EXIT_FAILURE);
	}

	run_over_blocks();

	return (EXIT_SUCCESS);
}

/******************************************************************************/

bool load_nc() {
	try {
		NcFile dataFile("../data/tos_O1_2001-2002.nc", NcFile::read);
		NcVar tosVar = dataFile.getVar("tos");
		if (tosVar.isNull()) {
			return false;
		}
		tosVar.getVar(data);

		for (int time = 0; time < NTIM; time++) {
			for (int lat = 0; lat < NLAT; lat++)
				for (int lon = 0; lon < NLON; lon++) {
					if (data[time][lat][lon] < 1000
							&& data[time][lat][lon] > -1000) {
//						printf("%f // tos(%d,%d,%d)\n", data[time][lat][lon],
//								time, lat, lon);
						data[time][lat][lon] = round(
								data[time][lat][lon] * 100000) / 100000.0f;
					} else {
//						printf("_ // tos(%d,%d,%d)\n", time, lat, lon);
						data[time][lat][lon] = 0;
					}
				}
		}
	} catch (NcException& e) {
		e.what();
		return false;
	}

	return true;
}

