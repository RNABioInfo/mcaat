# Use an official Ubuntu as a parent image
FROM ubuntu:latest
# Clean up unnecessary packages
# Install required packages
RUN apt-get update && apt-get -y upgrade
RUN apt-get install -y curl build-essential cmake make git pkg-config python3
RUN apt-get install -y libomp-dev zlib1g-dev libpthread-stubs0-dev 

# Set the working directory in the container
WORKDIR /mcaat

# Copy the project directory contents into the container at /mcaat
COPY . /mcaat
#CMD which
# Create build directory
RUN mkdir -p build

# Run CMake to configure the build
RUN cd build && cmake -DCMAKE_BUILD_TYPE=Release ..

# Build the project
RUN cd build && make VERBOSE=1

# Set the entry point for the container
ENTRYPOINT ["./build/mcaat"]
