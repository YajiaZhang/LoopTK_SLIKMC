################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
O_SRCS += \
../src/utils/camera/camera.o \
../src/utils/camera/clip.o \
../src/utils/camera/frustum.o \
../src/utils/camera/transform.o \
../src/utils/camera/viewport.o 

CPP_SRCS += \
../src/utils/camera/camera.cpp \
../src/utils/camera/clip.cpp \
../src/utils/camera/frustum.cpp \
../src/utils/camera/transform.cpp \
../src/utils/camera/viewport.cpp 

OBJS += \
./src/utils/camera/camera.o \
./src/utils/camera/clip.o \
./src/utils/camera/frustum.o \
./src/utils/camera/transform.o \
./src/utils/camera/viewport.o 

CPP_DEPS += \
./src/utils/camera/camera.d \
./src/utils/camera/clip.d \
./src/utils/camera/frustum.d \
./src/utils/camera/transform.d \
./src/utils/camera/viewport.d 


# Each subdirectory must supply rules for building sources it contributes
src/utils/camera/%.o: ../src/utils/camera/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cygwin C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


