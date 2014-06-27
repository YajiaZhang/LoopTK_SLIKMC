################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/utils/image/bmp.cpp \
../src/utils/image/gdi.cpp \
../src/utils/image/image.cpp \
../src/utils/image/import.cpp \
../src/utils/image/ppm.cpp \
../src/utils/image/textureops.cpp \
../src/utils/image/tga.cpp 

OBJS += \
./src/utils/image/bmp.o \
./src/utils/image/gdi.o \
./src/utils/image/image.o \
./src/utils/image/import.o \
./src/utils/image/ppm.o \
./src/utils/image/textureops.o \
./src/utils/image/tga.o 

CPP_DEPS += \
./src/utils/image/bmp.d \
./src/utils/image/gdi.d \
./src/utils/image/image.d \
./src/utils/image/import.d \
./src/utils/image/ppm.d \
./src/utils/image/textureops.d \
./src/utils/image/tga.d 


# Each subdirectory must supply rules for building sources it contributes
src/utils/image/%.o: ../src/utils/image/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cygwin C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


