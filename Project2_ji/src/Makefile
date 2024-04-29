# Compiler selection
CXX = g++

# Compiler flags
CXXFLAGS = -std=c++17 -Wall -Wextra

# Source files
SRCS = main.cpp Solver.cpp Mesh.cpp

# Header files
HDRS = Mesh.hpp Solver.hpp 

# Object files
OBJS = $(SRCS:.cpp=.o)

# Executable name
EXEC = Euler

# Targets
all: $(EXEC)

$(EXEC): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)

%.o: %.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(EXEC)
