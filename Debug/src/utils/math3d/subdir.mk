################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
O_SRCS += \
../src/utils/math3d/AABB2D.o \
../src/utils/math3d/AABB3D.o \
../src/utils/math3d/Circle3D.o \
../src/utils/math3d/Cylinder3D.o \
../src/utils/math3d/Line2D.o \
../src/utils/math3d/Line3D.o \
../src/utils/math3d/LinearlyDependent.o \
../src/utils/math3d/Plane3D.o \
../src/utils/math3d/Polygon2D.o \
../src/utils/math3d/Polyhedron3D.o \
../src/utils/math3d/Ray3D.o \
../src/utils/math3d/Segment2D.o \
../src/utils/math3d/Segment3D.o \
../src/utils/math3d/Triangle2D.o \
../src/utils/math3d/Triangle3D.o \
../src/utils/math3d/clip.o \
../src/utils/math3d/geometry2d.o \
../src/utils/math3d/geometry3d.o \
../src/utils/math3d/polar.o \
../src/utils/math3d/primitives.o \
../src/utils/math3d/quatinline.o \
../src/utils/math3d/random.o \
../src/utils/math3d/rotation.o 

CPP_SRCS += \
../src/utils/math3d/AABB2D.cpp \
../src/utils/math3d/AABB3D.cpp \
../src/utils/math3d/Circle3D.cpp \
../src/utils/math3d/Cylinder3D.cpp \
../src/utils/math3d/Line2D.cpp \
../src/utils/math3d/Line3D.cpp \
../src/utils/math3d/LinearlyDependent.cpp \
../src/utils/math3d/Plane3D.cpp \
../src/utils/math3d/Polygon2D.cpp \
../src/utils/math3d/Polyhedron3D.cpp \
../src/utils/math3d/Ray3D.cpp \
../src/utils/math3d/Segment2D.cpp \
../src/utils/math3d/Segment3D.cpp \
../src/utils/math3d/Triangle2D.cpp \
../src/utils/math3d/Triangle3D.cpp \
../src/utils/math3d/clip.cpp \
../src/utils/math3d/geometry2d.cpp \
../src/utils/math3d/geometry3d.cpp \
../src/utils/math3d/polar.cpp \
../src/utils/math3d/primitives.cpp \
../src/utils/math3d/quatinline.cpp \
../src/utils/math3d/random.cpp \
../src/utils/math3d/rotation.cpp 

OBJS += \
./src/utils/math3d/AABB2D.o \
./src/utils/math3d/AABB3D.o \
./src/utils/math3d/Circle3D.o \
./src/utils/math3d/Cylinder3D.o \
./src/utils/math3d/Line2D.o \
./src/utils/math3d/Line3D.o \
./src/utils/math3d/LinearlyDependent.o \
./src/utils/math3d/Plane3D.o \
./src/utils/math3d/Polygon2D.o \
./src/utils/math3d/Polyhedron3D.o \
./src/utils/math3d/Ray3D.o \
./src/utils/math3d/Segment2D.o \
./src/utils/math3d/Segment3D.o \
./src/utils/math3d/Triangle2D.o \
./src/utils/math3d/Triangle3D.o \
./src/utils/math3d/clip.o \
./src/utils/math3d/geometry2d.o \
./src/utils/math3d/geometry3d.o \
./src/utils/math3d/polar.o \
./src/utils/math3d/primitives.o \
./src/utils/math3d/quatinline.o \
./src/utils/math3d/random.o \
./src/utils/math3d/rotation.o 

CPP_DEPS += \
./src/utils/math3d/AABB2D.d \
./src/utils/math3d/AABB3D.d \
./src/utils/math3d/Circle3D.d \
./src/utils/math3d/Cylinder3D.d \
./src/utils/math3d/Line2D.d \
./src/utils/math3d/Line3D.d \
./src/utils/math3d/LinearlyDependent.d \
./src/utils/math3d/Plane3D.d \
./src/utils/math3d/Polygon2D.d \
./src/utils/math3d/Polyhedron3D.d \
./src/utils/math3d/Ray3D.d \
./src/utils/math3d/Segment2D.d \
./src/utils/math3d/Segment3D.d \
./src/utils/math3d/Triangle2D.d \
./src/utils/math3d/Triangle3D.d \
./src/utils/math3d/clip.d \
./src/utils/math3d/geometry2d.d \
./src/utils/math3d/geometry3d.d \
./src/utils/math3d/polar.d \
./src/utils/math3d/primitives.d \
./src/utils/math3d/quatinline.d \
./src/utils/math3d/random.d \
./src/utils/math3d/rotation.d 


# Each subdirectory must supply rules for building sources it contributes
src/utils/math3d/%.o: ../src/utils/math3d/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cygwin C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


