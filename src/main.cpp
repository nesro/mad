/*
 g++ -Ofast -march=native -Wall -Wextra -Wpedantic -std=c++11 -o main main.cpp -lnetcdf_c++ /home/n/Dropbox/netcdfcxx4/netcdf-cxx4-4.2.1/cxx4/.libs/libnetcdf_c++4.so.1 -lnetcdf -lbz2 bsdiff.c && time ./main
 */

#include <iostream>
#include <iomanip>

#include <functional>
#include <algorithm>
#include <tr1/functional>
#include <netcdf>
#include <string>
#include <bitset>
#include <tuple>
#include <cmath>
#include <cassert>
#include <limits>

#include "bzlib.h"
#include "bsdiff.h"

using namespace netCDF;
using namespace netCDF::exceptions;

static const int NLAT = 170;
static const int NLON = 180;
static const int NTIM = 24;

float data[NTIM][NLAT][NLON];

/* load data from tos_O1_2001-2002.nc to data */
bool load_nc();

static const int shl = (int) sizeof(size_t); /* simhash length */
static const int bs = 10; /* block size */
static const int binc = 4; /* block position increment */
static const int similarity_treshold = 10; /* how many bits can differ in simhash */
static const int search_window = 100;
static const float fp_eps = 0.001; /* test equality with: fabs(a,b) < fp_eps */
//static const int cout_precision = 10;
static const int close_together = 5;

void print_simhash(size_t simhash) {
	std::bitset<64> bs(simhash);
	std::cout << bs << std::endl;
}

/* stolen from boost */
void hash_combine(std::size_t& seed, const float& v) {
	std::tr1::hash<float> hasher;
	seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

size_t my_hash(float f) {
	std::tr1::hash<float> hasher;
	return (hasher(round(f * 1000)));
}

void cyclic_shift(size_t &s) {
	size_t first_bit = s << 63;
	s >>= 1;
	s |= first_bit;
}

short float2short(float f) {
	f *= 100;
	if (f > -32767 && f < 32767) {
		return ((short) f);
	} else {
		fprintf(stderr, "f too big :( %f\n", f);
		exit(1);
	}
}

float my_round(float f) {
	/* ncdump -p 9,17 -f c tos_O1_2001-2002.nc */
	return (roundf(f * 100.0f) / 100.0f);
}

typedef struct bsdiff_buf {
	unsigned char *buf;
	int size;
	int pos;
} bsdiff_buf_t;

int bsdiff_write(struct bsdiff_stream *stream, const void *buffer, int size) {
	bsdiff_buf_t *bsdiff_buf = (bsdiff_buf_t *) stream->opaque;

//	printf("bsdiff_write size = %d\n", size);

	if (bsdiff_buf->pos + size >= bsdiff_buf->size) {
		printf("not enough space\n");
		return (-1);
	}

	memcpy(bsdiff_buf->buf + bsdiff_buf->pos, buffer, size);
	bsdiff_buf->pos += size;

	return (0);
}

void block_bsdiff(int a1, int a2, int a3, int b1, int b2, int b3) {

	struct bsdiff_stream stream;

	const unsigned int buf_size = bs * bs * bs * sizeof(short);
	const unsigned int buf_size2 = 4096 * 8;

	unsigned char oldbuf[buf_size];
	short *soldbuf = (short *) oldbuf;
	int soldbuf_i = 0;

	unsigned char newbuf[buf_size];
	short *snewbuf = (short *) newbuf;
	int snewbuf_i = 0;

	unsigned char diff_buf[buf_size2];

	for (int d1 = 0; d1 < bs; d1++) {
		for (int d2 = 0; d2 < bs; d2++) {
			for (int d3 = 0; d3 < bs; d3++) {
				soldbuf[soldbuf_i++] = float2short(
						data[a1 + d1][a2 + d2][a3 + d3]);
				snewbuf[snewbuf_i++] = float2short(
						data[b1 + d1][b2 + d2][b3 + d3]);
			}
		}
	}

	stream.malloc = malloc;
	stream.free = free;
	stream.write = bsdiff_write;

	bsdiff_buf_t bsdiff_buf;
	bsdiff_buf.buf = diff_buf;
	bsdiff_buf.pos = 0;
	bsdiff_buf.size = buf_size2;

	stream.opaque = (void *) &bsdiff_buf;

	if (bsdiff(oldbuf, buf_size, newbuf, buf_size, &stream) == -1) {
		fprintf(stderr, "bsdiff failure\n");
		exit(1);
	}

	/* 1 - 9; 9 gives the best compression but uses the most runtime memory*/
	int blockSize = 9;
	/*1 - 4; 4 gives the most diagnostic info*/
	int verbosity = 0;
	/*30 is suggested; see docs for bzip2 for full info*/
	int workFactor = 30;

	char compr[buf_size2];
	unsigned int compr_size = buf_size2;

	int r = BZ2_bzBuffToBuffCompress(compr, &compr_size, (char *) diff_buf,
			bsdiff_buf.pos, blockSize, verbosity, workFactor);

	if (r != 0) {
		fprintf(stderr, "BZ2_bzBuffToBuffCompress r=%d\n", r);
		exit(1);
	}

	printf("bsdiff size=%d (r=%d)\n", compr_size, r);

}

unsigned int compress_block(int a1, int a2, int a3) {
	const unsigned int buf_size = bs * bs * bs * sizeof(short);
	char buf[buf_size];
	short *sbuf = (short *) buf;
	int sbuf_i = 0;

	char compr[buf_size];
	unsigned int compr_size = buf_size;

	for (int d1 = 0; d1 < bs; d1++) {
		for (int d2 = 0; d2 < bs; d2++) {
			for (int d3 = 0; d3 < bs; d3++) {
				if (0)
					printf("%d,", float2short(data[a1 + d1][a2 + d2][a3 + d3]));

				sbuf[sbuf_i++] = float2short(data[a1 + d1][a2 + d2][a3 + d3]);
			}
		}
	}

	/* 1 - 9; 9 gives the best compression but uses the most runtime memory*/
	int blockSize = 9;
	/*1 - 4; 4 gives the most diagnostic info*/
	int verbosity = 0;
	/*30 is suggested; see docs for bzip2 for full info*/
	int workFactor = 30;

	int r = BZ2_bzBuffToBuffCompress(compr, &compr_size, buf, buf_size,
			blockSize, verbosity, workFactor);
	printf("[%d][%d][%d] compr_size=%u r=%d\n", a1, a2, a3, compr_size, r);
	return (compr_size);
}

unsigned int compress_block_delta(int a1, int a2, int a3, int b1, int b2,
		int b3) {
	const unsigned int buf_size = bs * bs * bs * sizeof(short);
	char buf[buf_size];
	short *sbuf = (short *) buf;
	int sbuf_i = 0;

	char compr[buf_size];
	unsigned int compr_size = buf_size;

	for (int d1 = 0; d1 < bs; d1++) {
		for (int d2 = 0; d2 < bs; d2++) {
			for (int d3 = 0; d3 < bs; d3++) {

				if (0)
					printf("%d|",
							float2short(
									data[a1 + d1][a2 + d2][a3 + d3]
											- data[b1 + d1][b2 + d2][b3 + d3]));

				sbuf[sbuf_i++] = float2short(
						fabs(
								data[a1 + d1][a2 + d2][a3 + d3]
										- data[b1 + d1][b2 + d2][b3 + d3]));
			}
		}
	}

	/* 1 - 9; 9 gives the best compression but uses the most runtime memory*/
	int blockSize = 9;
	/*1 - 4; 4 gives the most diagnostic info*/
	int verbosity = 0;
	/*30 is suggested; see docs for bzip2 for full info*/
	int workFactor = 30;

	int r = BZ2_bzBuffToBuffCompress(compr, &compr_size, buf, buf_size,
			blockSize, verbosity, workFactor);
	printf("[%d][%d][%d] -> [%d][%d][%d] delta_compr_size=%u r=%d\n", a1, a2,
			a3, b1, b2, b3, compr_size, r);

	return (compr_size);
}

//size_t mdNNZsimhash(int d1, int d2, int d3) {
//	std::vector<int> mdsh_vec(64);
//	int nonzero = 0;
//
//	int nnz[2][2][2] = { 0 };
//
//	for (int i1 = 0; i1 < bs; i1++) {
//		for (int i2 = 0; i2 < bs; i2++) {
//			for (int i3 = 0; i3 < bs; i3++) {
//				float f = data[d1 + i1][d2 + i2][d3 + i3];
//
//				if (f == 0) {
//					continue;
//				}
//				nonzero++;
//				nnz[(i3 > 4)][(i2 > 4)][(i1 > 4)]++;
//
//				std::bitset<64> hash_bs(my_hash(f));
//
//				int index = 7
//						* (((i3 > 4) * 2 * 2) + ((i2 > 4) * 2) + (i1 > 4));
//
////				std::cout << i1 << "," << i2 << "," << i3 << "=" << f << ",ind="
////						<< index << std::endl;
//
//				for (int j = 0; j < 7; j++) {
//					if (hash_bs[j]) {
//						mdsh_vec[index + j]++;
//					} else {
//						mdsh_vec[index + j]--;
//					}
//				}
//
//			}
//		}
//	}
//
//	if (nonzero == 0) {
//		return (0);
//	}
//
//	std::bitset<64> mdsh(0);
//	for (int i = 0; i < 7 * 2 * 2 * 2; i++) {
//		if (mdsh_vec[i] >= 0) {
//			mdsh[i] = 1;
//		} else {
//			mdsh[i] = 0;
//		}
//	}
//
//	for (int i1 = 0; i1 < bs; i1++) {
//		for (int i2 = 0; i2 < bs; i2++) {
//			for (int i3 = 0; i3 < bs; i3++) {
//				if (nnz[i3][i2][i1] > 0) {
//					mdsh[7 * 2 * 2 * 2 + i3 * 4 + i2 * 2 + i1] = 1;
//				} else {
//					mdsh[7 * 2 * 2 * 2 + i3 * 4 + i2 * 2 + i1] = 0;
//				}
//			}
//		}
//	}
//
//	return ((size_t) mdsh.to_ulong());
//}

size_t mdsimhash(int d1, int d2, int d3) {
	std::vector<int> mdsh_vec(64);
	int nonzero = 0;

	for (int i1 = 0; i1 < bs; i1++) {
		for (int i2 = 0; i2 < bs; i2++) {
			for (int i3 = 0; i3 < bs; i3++) {
				float f = data[d1 + i1][d2 + i2][d3 + i3];

				if (f == 0) {
					continue;
				}
				nonzero++;

				std::bitset<64> hash_bs(my_hash(f));

				int index = 8
						* (((i3 > 4) * 2 * 2) + ((i2 > 4) * 2) + (i1 > 4));

//				std::cout << i1 << "," << i2 << "," << i3 << "=" << f << ",ind="
//						<< index << std::endl;

				for (int j = 0; j < 8; j++) {
					if (hash_bs[j]) {
						mdsh_vec[index + j]++;
					} else {
						mdsh_vec[index + j]--;
					}
				}

			}
		}
	}

	if (nonzero < (bs * bs * bs / 2)) {
		return (0);
	}

	std::bitset<64> mdsh(0);
	for (int i = 0; i < 64; i++) {
		if (mdsh_vec[i] >= 0) {
			mdsh[i] = 1;
		} else {
			mdsh[i] = 0;
		}
	}

	return ((size_t) mdsh.to_ulong());
}

/******************************************************************************/

class Point {
public:
	int x;
	int y;
	int z;
	float v;

	Point(int _x, int _y, int _z, float _v) :
			x(_x), y(_y), z(_z), v(_v) {
	}
};

class Block {
public:
	int tim;
	int lat;
	int lon;
	std::size_t simhash;
	std::vector<Block> similar;

	Block(int t, int la, int lo, std::size_t sh) :
			tim(t), lat(la), lon(lo), simhash(sh) {
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

void print_block_data(Block &b1);

/******************************************************************************/

/* Give block of 8x8x8 values.
 * Return a 64bit simhash. It represents a 2x2x2 cube with 8bit values.
 [00] [01] [02] [03] [04] [05] [06] [07] - (0,0,0)
 [08] [09] [10] [11] [12] [13] [14] [15] - (0,0,1)
 [16] [17] [18] [19] [20] [21] [22] [23] - (0,1,0)
 [24] [25] [26] [27] [28] [29] [30] [31] - (0,1,1)
 [32] [33] [34] [35] [36] [37] [38] [39] - (1,0,0)
 [40] [41] [42] [43] [44] [45] [46] [47] - (1,0,1)
 [48] [49] [50] [51] [52] [53] [54] [55] - (1,1,0)
 [56] [57] [58] [59] [60] [61] [62] [63] - (1,1,1)
 */

/******************************************************************************/

class SimilarBlocks {
public:
	Block b1;
	Block b2;
	int zero_same;
	int same;
	float diff;
	float diff_avg;
	int xordiff;
	std::vector<Point> diffpoints;

	SimilarBlocks(Block &aa, Block &bb) :
			b1(aa), b2(bb) {

		zero_same = 0;
		same = 0;
		diff = 0.0f;
		xordiff = 0;

		for (int iti = 0; iti < bs; iti++) {
			for (int ilat = 0; ilat < bs; ilat++) {
				for (int ilon = 0; ilon < bs; ilon++) {
					float fb1 = data[b1.tim + iti][b1.lat + ilat][b1.lon + ilon];
					float fb2 = data[b2.tim + iti][b2.lat + ilat][b2.lon + ilon];

					if (fabs(fb1 - fb2) < fp_eps) {
						if (fb1 == 0) {
							zero_same++;
							continue;
						}
						same++;
						continue;
					}

					diffpoints.push_back(
							Point(b1.tim + iti, b1.lat + ilat, b1.lon + ilon,
									my_round(fb1 - fb2)));

					diff += fabs(my_round(fb1 - fb2));
				}
			}
		}

		diff_avg = diff / diffpoints.size();

		std::bitset<64> simxor(b1.simhash ^ b2.simhash);
		xordiff = simxor.count();

	}

	bool operator <(const SimilarBlocks& other) const {
		return (this->diff_avg < other.diff_avg);
	}

	void print() {
		std::cout << "b1[" << b1.lat << "][" << b1.lon << "][" << b1.tim
				<< "], ";
		std::cout << "b2[" << b2.lat << "][" << b2.lon << "][" << b2.tim
				<< "], ";

		std::cout << "zero_same=" << zero_same << ", same=" << same << ", diff="
				<< diff << ", diff_avg=" << diff_avg << ", xor=" << xordiff
				<< std::endl;
		print_simhash(b1.simhash);
		print_simhash(b2.simhash);

		std::cout << "print b1" << std::endl;
		print_block_data(b1);
		std::cout << std::endl;

		std::cout << "print b2" << std::endl;
		print_block_data(b2);
		std::cout << std::endl;

//		std::cout.precision(std::numeric_limits<double>::digits10);
//		std::cout.setf(std::ios::fixed, std::ios::floatfield);

		std::cout << "print diff" << std::endl;
		for (auto &d : diffpoints) {
			std::cout << std::setprecision(7) << d.v << ", ";
		}
		std::cout << std::endl;
		std::cout << std::endl;
	}
};

/******************************************************************************/

class BlockDelta {
public:
	Block &from;
	Block &to;

	BlockDelta(Block &_from, Block &_to) :
			from(_from), to(_to) {

	}

};

/******************************************************************************/

void print_block_data(Block &b1) {
	for (int iti = 0; iti < bs; iti++) {
		for (int ilat = 0; ilat < bs; ilat++) {
			for (int ilon = 0; ilon < bs; ilon++) {
//				std::cout.setf(std::ios::fixed, std::ios::floatfield);
				std::cout.precision(12);
				std::cout << data[b1.tim + iti][b1.lat + ilat][b1.lon + ilon]
						<< ", ";
			}
		}
	}
	std::cout << std::endl;
}

/******************************************************************************/

std::vector<SimilarBlocks> find_similar(std::vector<Block> simhashes) {
	std::vector<SimilarBlocks> similarBlocks;

	for (int s = 0; s < 64; s++) {
		std::sort(simhashes.begin(), simhashes.end());

		for (int i = 0; i < (int) simhashes.size() - search_window; i++) {
			for (int sw = 0; sw < search_window; sw++) {
				std::bitset<64> simxor(
						simhashes[i].simhash ^ simhashes[i + sw].simhash);

				if (simxor.count() < similarity_treshold) {
					/* skip blocks too close together */
					if ((abs(simhashes[i].tim - simhashes[i + sw].tim)
							< close_together)
							|| (abs(simhashes[i].lat - simhashes[i + sw].lat)
									< close_together)
							|| (abs(simhashes[i].lon - simhashes[i + sw].lon)
									< close_together)) {
						continue;
					}

					/* are they already similar? */
					if (std::find(simhashes[i].similar.begin(),
							simhashes[i].similar.end(), simhashes[i + sw])
							!= simhashes[i].similar.end()) {
						continue;
					}

					simhashes[i].similar.push_back(simhashes[i + sw]);
					simhashes[i + sw].similar.push_back(simhashes[i]);

					similarBlocks.push_back(
							SimilarBlocks(simhashes[i], simhashes[i + sw]));
				}
			}
		}

		for (int i = 0; i < (int) simhashes.size(); i++) {
			cyclic_shift(simhashes[i].simhash);
		}
	}

	/* sort them by differcence */
	std::cout << "size: " << similarBlocks.size() << ", now sort" << std::endl;
	std::sort(similarBlocks.begin(), similarBlocks.end());

	return (similarBlocks);
}

void run_over_blocks2() {
	std::vector<Block> simhashes;

	for (int tim = 0; tim < NTIM - bs; tim += binc) {
		for (int lat = 0; lat < NLAT - bs; lat += binc) {
			for (int lon = 0; lon < NLON - bs; lon += binc) {
				size_t sh = mdsimhash(tim, lat, lon);
				if (sh > 0) {
					simhashes.push_back(Block(tim, lat, lon, sh));
				}
			}
		}
	}
	std::cout << "simhashes.sie() = " << simhashes.size() << std::endl;

	std::cout << "find similar blocks" << std::endl;
	std::vector<SimilarBlocks> similarBlocks = find_similar(simhashes);

	std::cout << "found similar: " << similarBlocks.size() << std::endl;

	if (similarBlocks.size() == 0) {
		std::cout << "no similar. end" << std::endl;
		return;
	}

	for (int i = 0; i < 5; i++) {
		std::cout << "i=" << i << std::endl;
		similarBlocks[i].print();
	}
}

//void run_over_blocks() {
//
//	std::tr1::hash<int> hasher;
//	std::vector<Block> simhashes;
//	std::vector<std::pair<Block, Block>> similar;
//	std::vector<SimilarBlocks> similarBlocks;
//
//	for (int tim = 0; tim < NTIM - bs; tim += binc) {
//		for (int lat = 0; lat < NLAT - bs; lat += binc) {
//			for (int lon = 0; lon < NLON - bs; lon += binc) {
//				std::vector<int> simhash_vec(64);
//				int nonzero = 0;
//
//				for (int iti = 0; iti < bs; iti++) {
//					for (int ilat = 0; ilat < bs; ilat++) {
//						for (int ilon = 0; ilon < bs; ilon++) {
//
//							std::size_t seed = 0;
//							float f = data[tim + iti][lat + ilat][lon + ilon];
//
//							if (f == 0) {
//								continue;
//							}
//							nonzero++;
//
//							// https://stackoverflow.com/questions/19966041/getting-too-many-collisions-with-hash-combine
//							hash_combine(seed, hasher(iti) * 2654435761);
//							hash_combine(seed,
//									hasher(ilat * 1000) * 2654435761);
//							hash_combine(seed,
//									hasher(ilon * 1000000) * 2654435761);
////							hash_combine(seed,
////									hasher(hasher(tim + lat * 1000 + lon * 100000))
////											* 2654435761);
//							hash_combine(seed, my_hash(f) * 2654435761);
//
//							std::bitset<64> seed_bitset(seed);
////					std::cout << "seed_bitset: " << seed_bitset << std::endl;
//							for (int j = 0; j < 64; j++) {
//								if (seed_bitset[j]) {
//									simhash_vec[j]++;
//								} else {
//									simhash_vec[j]--;
//								}
//							}
//
////					std::cout << "simhash_vec: ";
////					for (auto i = simhash_vec.begin(); i != simhash_vec.end(); ++i)
////					    std::cout << *i << ", ";
////					std::cout << std::endl;
//
//						}
//					}
//				} /* for i 0..bs */
//
//				/* if there was no nonzero value, don't save its simhash */
//				if (nonzero == 0) {
//					continue;
//				}
//
//				std::bitset<shl> simhash_bitset(0);
//				for (int i = 0; i < 64; i++) {
//					if (simhash_vec[i] >= 0) {
//						simhash_bitset[i] = 1;
//					} else {
//						simhash_bitset[i] = 0;
//					}
//				}
//
//				std::size_t simhash = simhash_bitset.to_ulong();
//				simhashes.push_back(Block(tim, lat, lon, simhash));
//			}
//		}
//	}
//
//	std::cout << "now sort" << std::endl;
//
//	for (int s = 0; s < 64; s++) {
//		std::sort(simhashes.begin(), simhashes.end());
//
//		for (int i = 0; i < (int) simhashes.size() - search_window; i++) {
//			for (int sw = 0; sw < search_window; sw++) {
//				std::bitset<64> simxor(
//						simhashes[i].simhash ^ simhashes[i + sw].simhash);
//
//				if (simxor.count() < similarity_treshold) {
//
//					/* skip blocks too close together */
//					/* TODO: make 5 b1 constant */
//					if ((abs(simhashes[i].tim - simhashes[i + sw].tim) < 5)
//							|| (abs(simhashes[i].lat - simhashes[i + sw].lat)
//									< 5)
//							|| (abs(simhashes[i].lon - simhashes[i + sw].lon)
//									< 5)) {
//						continue;
//					}
//
//					/* are they already similar? */
//					if (std::find(simhashes[i].similar.begin(),
//							simhashes[i].similar.end(), simhashes[i + sw])
//							!= simhashes[i].similar.end()) {
//						continue;
//					}
//
////					std::cout << "found almost same simhash" << std::endl;
////
////					std::cout << "  first: ";
////					std::cout << "tim=" << simhashes[i].tim << ", ";
////					std::cout << "lon=" << simhashes[i].lon << ", ";
////					std::cout << "lat=" << simhashes[i].lat;
////					std::cout << std::endl;
////					std::cout << "  ";
////					print_simhash(simhashes[i].simhash);
////
////					std::cout << "  second: ";
////					std::cout << "tim=" << simhashes[i + sw].tim << ", ";
////					std::cout << "lon=" << simhashes[i + sw].lon << ", ";
////					std::cout << "lat=" << simhashes[i + sw].lat;
////					std::cout << std::endl;
////					std::cout << "  ";
////						simhash(simhashes[i + sw].simhash);
////
////					std::cout << "  ";
////					block_diff(simhashes[i], simhashes[i + sw]);
////					std::cout << std::endl;
//
//					simhashes[i].similar.push_back(simhashes[i + sw]);
//					simhashes[i + sw].similar.push_back(simhashes[i]);
//
//					similar.push_back(
//							std::make_pair(simhashes[i], simhashes[i + sw]));
//
//					similarBlocks.push_back(
//							SimilarBlocks(simhashes[i], simhashes[i + sw]));
//				}
//			}
//		}
//
//		for (int i = 0; i < (int) simhashes.size(); i++) {
//			cyclic_shift(simhashes[i].simhash);
//		}
//	}
//
//	std::cout << "sorting similars" << std::endl;
//
//	std::sort(similarBlocks.begin(), similarBlocks.end());
//
//	std::cout << "found similar: " << similarBlocks.size() << std::endl;
//
//	for (int i = 0; i < 1; i++) {
//		std::cout << "i=" << i << std::endl;
//		similarBlocks[i].print();
//	}
//
//}

/******************************************************************************/

void print_block(int a1, int a2, int a3) {
	printf("{[%d][%d][%d]:\n", a1, a2, a3);
	for (int d1 = 0; d1 < bs; d1++) {
		printf("[");
		for (int d2 = 0; d2 < bs; d2++) {
			printf("(");
			for (int d3 = 0; d3 < bs; d3++) {
				float f = data[a1 + d1][a2 + d2][a3 + d3];
				if (f > 0) {
					printf("%.2f,", f);
				} else {
					printf(",");
				}
			}
			printf(")");
		}
		printf("]\n");
	}
	printf("}\n");
}

/* return sum */
void print_block_diff(int a1, int a2, int a3, int b1, int b2, int b3,
		bool print, double &diff, int &diff_cnt, int &non_zero) {
	if (print) {
		printf("{[%d][%d][%d]delta[%d][%d][%d]:\n", a1, a2, a3, b1, b2, b3);
	}
	for (int d1 = 0; d1 < bs; d1++) {
		if (print) {
			printf("[");
		}
		for (int d2 = 0; d2 < bs; d2++) {
			if (print) {
				printf("(");
			}
			for (int d3 = 0; d3 < bs; d3++) {
				float f1 = data[a1 + d1][a2 + d2][a3 + d3];
				float f2 = data[b1 + d1][b2 + d2][b3 + d3];

				if (f1 > 0 || f2 > 0) {
					non_zero++;
				}

				float f3 = f1 - f2;

				if (fabs(f3) > fp_eps) {
					diff_cnt++;
					diff += fabs(f3);
				}

				if (print) {
					if (fabs(f3) > fp_eps) {
						printf("%.2f,", f3);
					} else {
						printf(",");
					}
				}
			}
			if (print) {
				printf(")");
			}
		}
		if (print) {
			printf("]\n");
		}
	}
	if (print) {
		printf(":sum=%f}\n", diff);
	}
}

//class Sim {
//public:
//	int size1;
//	int size2;
//	int size_subtract_diff;
//	int size_bsdiff_diff;
//	int better;
//
//	Sim(int s1, int s2, int ss, int sb) :
//			size1(s1), size2(s2), size_subtract_diff(ss), size_bsdiff_diff(sb) {
//
//		int ssx = size2 - size_subtract_diff;
//		int ssy = size2 - size_bsdiff_diff;
//
//		better = std::max(ssx, ssy);
//
//	}
//
//	bool operator <(const Sim& other) const {
//			return (this->size2 < other.diff_avg);
//		}
//};

void brute_force() {
	int similar = 0;
	for (int tim = 0; tim < NTIM - bs; tim += binc) {
		printf("tim=%d\n", tim);
		for (int lat = 0; lat < NLAT - bs; lat += binc) {
			for (int lon = 0; lon < NLON - bs; lon += binc) {

				for (int tim2 = tim + close_together; tim2 < NTIM - bs; tim2 +=
						binc) {
					for (int lat2 = lat + close_together; lat2 < NLAT - bs;
							lat2 += binc) {
						for (int lon2 = lon + binc + close_together;
								lon2 < NLON - bs; lon2 += binc) {

							double diff = 0;
							int diff_cnt = 0;
							int non_zero = 0;

							print_block_diff(tim, lat, lon, tim2, lat2, lon2,
									false, diff, diff_cnt, non_zero);
							double diff_avg = diff / (float) diff_cnt;

//							printf("diff=%f diff_cnt=%d diff_avg=%f\n", diff,
//									diff_cnt, diff_avg);

							if (non_zero > 0 && fabs(diff_avg) < 1) {
								printf(
										"diff=%03f,\tdiff_avg=%f\t[%03d][%03d][%03d] [%03d][%03d][%03d]\n",
										diff, diff_avg, tim, lat, lon, tim2,
										lat2, lon2);
								if (1) {
									if (1) {
										print_block(tim, lat, lon);
										print_block(tim2, lat2, lon2);
										print_block_diff(tim, lat, lon, tim2,
												lat2, lon2, true, diff,
												diff_cnt, non_zero);
									}

									compress_block(tim, lat, lon);
									compress_block(tim2, lat2, lon2);
									compress_block_delta(tim, lat, lon, tim2,
											lat2, lon2);
									block_bsdiff(tim2, lat2, lon2, tim, lat,
											lon);
									block_bsdiff(tim, lat, lon, tim2, lat2,
											lon2);

									printf("----\n\n");
								}
								similar++;
							}
						}
					}
				}

			}

		}
	}
	printf("total similar=%d\n", similar);
}

/******************************************************************************/

/******************************************************************************/

int main() {

//	std::cout << "h1: " << my_hash(1.0001f) << std::endl;
//	std::cout << "h1: " << my_hash(1.0002f) << std::endl;
//	return (EXIT_SUCCESS);

	if (!load_nc()) {
		return (EXIT_FAILURE);
	}

//	compress_block(4, 152, 843);
//	compress_block(9, 157, 93);
//	compress_block_delta(9, 157, 93, 4, 152, 843);
//	compress_block_delta(4, 152, 843, 9, 157, 93);

//	printf("bsdiff\n");
//	block_bsdiff(9, 157, 93, 4, 152, 843);
//	block_bsdiff(4, 152, 843, 9, 157, 93);

	brute_force();

//	run_over_blocks2();

//	std::cout << mdsimhash(15, 25, 53) << std::endl;

//	run_over_blocks();

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

		int nnz = 0;
		int zero = 0;

		for (int time = 0; time < NTIM; time++) {
			for (int lat = 0; lat < NLAT; lat++)
				for (int lon = 0; lon < NLON; lon++) {
					if (data[time][lat][lon] < 1000
							&& data[time][lat][lon] > -1000) {
//						printf("%f // tos(%d,%d,%d)\n", data[time][lat][lon],
//								time, lat, lon);
						data[time][lat][lon] = my_round(data[time][lat][lon]);
						nnz++;
					} else {
//						printf("_ // tos(%d,%d,%d)\n", time, lat, lon);
						data[time][lat][lon] = 0;
						zero++;
					}
				}
		}
		printf("nnz=%d, zero=%d\n", nnz, zero);
	} catch (NcException& e) {
		e.what();
		return false;
	}

//	std::cout.precision(12);
//	std::cout << data[0][78][77] << std::endl;

	return true;
}

