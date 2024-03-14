#include "utils.hpp"

int utils::stoe(std::string& s) {
    int id = -1;

    if(s == "Mono16")
        id = Mono16;
    if(s == "Mono12")
        id = Mono12;
    if(s == "Mono12Packed")
        id = Mono12Packed;
    if(s == "Mono32")
        id = Mono32;

    return id;
}

fftwf_complex* utils::allocate_complex_float_array(size_t N) {
#ifdef __AVX512F__
	return (fftwf_complex*) _mm_malloc(2*N*sizeof(float), 64);
#else
	return fftwf_alloc_complex(N);
#endif
}

float* utils::allocate_float_array(size_t N)  {
#ifdef __AVX512F__
	return (float*) _mm_malloc(N*sizeof(float), 64);
#else
	return fftwf_alloc_real(N);
#endif
}

template <typename ImType>
bool utils::libTIFFWriter_writeImage(
		TIFF* tif, ImType* img, int width, int height) {

	uint16_t frame_width = width;
	uint16_t frame_height = height;
	uint32_t rowsperstrip = (uint32_t)-1;

	TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, frame_width);
	TIFFSetField(tif, TIFFTAG_IMAGELENGTH, frame_height);
	TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
	TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
	TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
	TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, sizeof(ImType)*8);
	TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 1);

	rowsperstrip = TIFFDefaultStripSize(tif, rowsperstrip);
	TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, rowsperstrip);
	TIFFSetField(tif, TIFFTAG_FILLORDER, FILLORDER_MSB2LSB);
	TIFFSetField(tif,TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);
	TIFFSetField(tif,TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);

	// Data is broken up into strips where each strip contains rowsperstrip complete rows of
	// data. The last strip is possibly smaller, so it is NOT padded with dummy data.
	for(unsigned int row = 0; row<frame_height; row+=rowsperstrip) {
		// compute rows in this strip:
		uint32_t nrow = rowsperstrip>frame_height-row ? frame_height-row : rowsperstrip;
		tstrip_t strip = TIFFComputeStrip(tif, row, 0);
		if(TIFFWriteEncodedStrip(
					tif, strip, &img[row*frame_width], sizeof(ImType)*nrow*frame_width)<0)
			return false;
	}
	TIFFWriteDirectory(tif);
	return true;
}

template bool utils::libTIFFWriter_writeImage<float>(
		TIFF* tif, float* img, int width, int height);