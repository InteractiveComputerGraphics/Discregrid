#pragma once

class BmpReaderWriter
{
public:
	static bool isBigEndian();
	static unsigned short endianSwap(unsigned short nValue);
	static unsigned int endianSwap(unsigned int i);

	// -------------------------------------------------------------------

#pragma pack(1)

	struct BMPHEADER {
		unsigned short		Type;
		unsigned int		Size;
		unsigned short		Reserved1;
		unsigned short		Reserved2;
		unsigned int		OffBits;
	};

	// Only Win3.0 BMPINFO (see later for OS/2)
	struct BMPINFO {
		unsigned int		Size;
		unsigned int		Width;
		unsigned int		Height;
		unsigned short		Planes;
		unsigned short		BitCount;
		unsigned int		Compression;
		unsigned int		SizeImage;
		unsigned int		XPelsPerMeter;
		unsigned int		YPelsPerMeter;
		unsigned int		ClrUsed;
		unsigned int		ClrImportant;
	};

#pragma pack()

	// Data is persists until the class is destructed.
	static bool loadFile(const char *filename, unsigned int &width, unsigned int &height, unsigned char *&data);
	static bool saveFile(const char *filename, int width, int height, unsigned char *data);
};
