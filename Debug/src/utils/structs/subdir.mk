################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/utils/structs/set.cpp 

OBJS += \
./src/utils/structs/set.o 

CPP_DEPS += \
./src/utils/structs/set.d 


# Each subdirectory must supply rules for building sources it contributes
src/utils/structs/%.o: ../src/utils/structs/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cygwin C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


