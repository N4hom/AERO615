# Compiler selection
CXX = g++

# Compiler flags
CXXFLAGS = -std=c++11 -Wall -Wextra

# Source files
SRCS = main.cpp Problem.cpp 

# Header files
HDRS = Problem.hpp Matrix.hpp Vector.hpp

# Object files
OBJS = $(SRCS:.cpp=.o)

# Executable name
EXEC = gridGen

# Targets
all: $(EXEC)

$(EXEC): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)

%.o: %.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(EXEC)
