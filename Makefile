# Load CUDA using the following command
# module load cuda
#
#CC = nvcc
#CFLAGS = -O3 -arch=compute_37 -code=sm_37
#NVCCFLAGS = -O3 -arch=compute_37 -code=sm_37
#LIBS = 

CC = nvcc
CFLAGS = -O3 -arch=compute_61 -code=sm_61
NVCCFLAGS = -O3 -arch=compute_61 -code=sm_61
LIBS = 

TARGETS = serial gpu autograder

all:	$(TARGETS)

serial: serial.o common.o core.o
	$(CC) -o $@ $(LIBS) serial.o common.o core.o
gpu: gpu.o common.o
	$(CC) -o $@ $(NVCCLIBS) gpu.o common.o core.o
autograder: autograder.o common.o
	$(CC) -o $@ $(LIBS) autograder.o common.o

serial.o: serial.cu common.h core.h
	$(CC) -c $(CFLAGS) serial.cu
autograder.o: autograder.cu common.h
	$(CC) -c $(CFLAGS) autograder.cu
gpu.o: gpu.cu common.h core.h
	$(CC) -c $(NVCCFLAGS) gpu.cu
common.o: common.cu common.h
	$(CC) -c $(CFLAGS) common.cu
core.o: core.cu core.h
	$(CC) -c $(CFLAGS) core.cu

clean:
	rm -f *.o $(TARGETS) *.stdout
