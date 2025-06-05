# Makefile

CXX = g++
CXXFLAGS = -std=c++17 -O2


INCLUDES = -I/opt/homebrew/opt/hdf5/include
LIBS = -L/opt/homebrew/opt/hdf5/lib -lhdf5_cpp -lhdf5  # <-- if you downloaded HighFive headers, otherwise skip
#LIBS = -lhdf5_cpp -lhdf5

TARGET = alphadynamics_128_0pt95
SRCS = alphadynamics_2.cpp

all: $(TARGET)

$(TARGET): $(SRCS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(LIBS) -o $(TARGET) $(SRCS)

clean:
	rm -f $(TARGET)


