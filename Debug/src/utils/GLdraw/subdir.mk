################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
O_SRCS += \
../src/utils/GLdraw/GLColor.o \
../src/utils/GLdraw/GLError.o \
../src/utils/GLdraw/GLFog.o \
../src/utils/GLdraw/GLLight.o \
../src/utils/GLdraw/GLMaterial.o \
../src/utils/GLdraw/GLScreenshot.o \
../src/utils/GLdraw/GLTextureObject.o \
../src/utils/GLdraw/GLUINavigationProgram.o \
../src/utils/GLdraw/GLUIProgram.o \
../src/utils/GLdraw/GLUTNavigationProgram.o \
../src/utils/GLdraw/GLUTProgram.o \
../src/utils/GLdraw/GLView.o \
../src/utils/GLdraw/drawextra.o 

CPP_SRCS += \
../src/utils/GLdraw/GLColor.cpp \
../src/utils/GLdraw/GLError.cpp \
../src/utils/GLdraw/GLFog.cpp \
../src/utils/GLdraw/GLLight.cpp \
../src/utils/GLdraw/GLMaterial.cpp \
../src/utils/GLdraw/GLScreenshot.cpp \
../src/utils/GLdraw/GLTextureObject.cpp \
../src/utils/GLdraw/GLUINavigationProgram.cpp \
../src/utils/GLdraw/GLUIProgram.cpp \
../src/utils/GLdraw/GLUTNavigationProgram.cpp \
../src/utils/GLdraw/GLUTProgram.cpp \
../src/utils/GLdraw/GLView.cpp \
../src/utils/GLdraw/drawextra.cpp 

OBJS += \
./src/utils/GLdraw/GLColor.o \
./src/utils/GLdraw/GLError.o \
./src/utils/GLdraw/GLFog.o \
./src/utils/GLdraw/GLLight.o \
./src/utils/GLdraw/GLMaterial.o \
./src/utils/GLdraw/GLScreenshot.o \
./src/utils/GLdraw/GLTextureObject.o \
./src/utils/GLdraw/GLUINavigationProgram.o \
./src/utils/GLdraw/GLUIProgram.o \
./src/utils/GLdraw/GLUTNavigationProgram.o \
./src/utils/GLdraw/GLUTProgram.o \
./src/utils/GLdraw/GLView.o \
./src/utils/GLdraw/drawextra.o 

CPP_DEPS += \
./src/utils/GLdraw/GLColor.d \
./src/utils/GLdraw/GLError.d \
./src/utils/GLdraw/GLFog.d \
./src/utils/GLdraw/GLLight.d \
./src/utils/GLdraw/GLMaterial.d \
./src/utils/GLdraw/GLScreenshot.d \
./src/utils/GLdraw/GLTextureObject.d \
./src/utils/GLdraw/GLUINavigationProgram.d \
./src/utils/GLdraw/GLUIProgram.d \
./src/utils/GLdraw/GLUTNavigationProgram.d \
./src/utils/GLdraw/GLUTProgram.d \
./src/utils/GLdraw/GLView.d \
./src/utils/GLdraw/drawextra.d 


# Each subdirectory must supply rules for building sources it contributes
src/utils/GLdraw/%.o: ../src/utils/GLdraw/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cygwin C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


