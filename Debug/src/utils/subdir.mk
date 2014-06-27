################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
O_SRCS += \
../src/utils/Timer.o \
../src/utils/myfile.o 

CPP_SRCS += \
../src/utils/Timer.cpp \
../src/utils/myfile.cpp 

OBJS += \
./src/utils/Timer.o \
./src/utils/myfile.o 

CPP_DEPS += \
./src/utils/Timer.d \
./src/utils/myfile.d 


# Each subdirectory must supply rules for building sources it contributes
src/utils/%.o: ../src/utils/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cygwin C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


