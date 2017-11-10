#include <cstdio>
#include <cstdlib>
#include "bmp_file.hpp"

// Compression Type
#define BI_RGB      0L
#define BI_RLE8     1L
#define BI_RLE4     2L


// -------------------------------------------------------------------
bool BmpReaderWriter::loadFile(const char *filename, unsigned int &width, unsigned int &height, unsigned char *&data)
{
	if (data) 
	{		
		delete [] data;
		data = NULL;
	}
	width = 0;
	height = 0;

	FILE *f = fopen(filename, "rb");
	if (!f) 
		return false;

	size_t num;
	BMPHEADER header;
	num = fread(&header, sizeof(BMPHEADER), 1, f);
	if(isBigEndian()) header.Type = endianSwap(header.Type);
	if (num != 1) { fclose(f); return false; }
	if (header.Type != 'MB') { fclose(f); return false; }

	BMPINFO info;
	num = fread(&info, sizeof(BMPINFO), 1, f);
	if (num != 1) { fclose(f); return false; }
	if(isBigEndian()) info.Size = endianSwap(info.Size);
	if(isBigEndian()) info.BitCount = endianSwap(info.BitCount);
	if(isBigEndian()) info.Compression = endianSwap(info.Compression);
	if(isBigEndian()) info.Width = endianSwap(info.Width);
	if(isBigEndian()) info.Height = endianSwap(info.Height);

	if (info.Size != sizeof(BMPINFO)) { fclose(f); return false; }
	if (info.BitCount != 24) { fclose(f); return false; }
	if (info.Compression != BI_RGB) { fclose(f); return false; }

	width = info.Width;
	height = info.Height;
	data = new unsigned char[width * height * 3];

	int lineLen = (((info.Width * (info.BitCount>>3)) + 3)>>2)<<2;
	unsigned char *line = new unsigned char[lineLen];

	for(int i = info.Height-1; i >= 0; i--) {
		num = fread(line, lineLen, 1, f);
		if (num != 1) { fclose(f); return false; }
		unsigned char *src = line;
		unsigned char *dest = data + i*info.Width*3;
		for(unsigned int j = 0; j < info.Width; j++) {
			unsigned char r,g,b;
			b = *src++; g = *src++; r = *src++;
			*dest++ = r; *dest++ = g; *dest++ = b;
		}
	}

	delete [] line;
	fclose(f);

	return true;
}

// -------------------------------------------------------------------
bool BmpReaderWriter::saveFile(const char *filename, int width, int height, unsigned char *data)
{
	FILE *f = fopen(filename, "wb");
	if (!f) return false;

	// todo : works on pcs only, swap correctly if big endian 
	BMPHEADER header;
	header.Type = 'MB';
	header.Size = sizeof(BMPINFO);
	header.Reserved1 = 0;
	header.Reserved2 = 0;
	header.OffBits = sizeof(BMPHEADER) + sizeof(BMPINFO);
	fwrite(&header, sizeof(BMPHEADER), 1, f);

	BMPINFO info;
	info.Size = sizeof(BMPINFO);
	info.Width = width;
	info.Height = height;
	info.Planes = 1;
	info.BitCount = 24;
	info.Compression = BI_RGB;
	info.XPelsPerMeter = 4000;
	info.YPelsPerMeter = 4000;
	info.ClrUsed = 0;
	info.ClrImportant = 0;
	fwrite(&info, sizeof(info), 1, f);

	// padded to multiple of 4
	int lineLen = (((info.Width * (info.BitCount>>3)) + 3)>>2)<<2;
	info.SizeImage = lineLen * height;

	unsigned char *line = new unsigned char[lineLen];

	for(int i = 0; i < height; i++) 
	{
		unsigned char *src = data + i*width*3;
		unsigned char *dest = line;
		for(int j = 0; j < width; j++) 
		{
			unsigned char r,g,b;
			r = *src++; g = *src++; b = *src++;
			*dest++ = b; *dest++ = g; *dest++ = r;
		}
		for (int j = 3*width; j < lineLen; j++)
			*dest++ = 0;
		fwrite(line, lineLen, 1, f);
	}

	delete [] line;
	fclose(f);

	return true;
}

bool BmpReaderWriter::isBigEndian()
{
	int i = 1; return *((char*)&i) == 0;
}

unsigned short BmpReaderWriter::endianSwap(unsigned short nValue)
{
	return (((nValue >> 8)) | (nValue << 8));
}

unsigned int BmpReaderWriter::endianSwap(unsigned int i)
{
	unsigned char b1, b2, b3, b4;

	b1 = i & 255;
	b2 = (i >> 8) & 255;
	b3 = (i >> 16) & 255;
	b4 = (i >> 24) & 255;

	return ((unsigned int)b1 << 24) + ((unsigned int)b2 << 16) + ((unsigned int)b3 << 8) + b4;
}

