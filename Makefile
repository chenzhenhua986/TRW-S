all : stereo

stereo: CImg.h loopy_BP.cpp
	g++ -Dcimg_display=0 loopy_BP.cpp -o loopy_BP -I. -O3
clean:
	rm loopy_BP
