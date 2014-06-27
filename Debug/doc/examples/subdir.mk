################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
O_SRCS += \
../doc/examples/main.o 

CC_SRCS += \
../doc/examples/main.cc 

OBJS += \
./doc/examples/main.o 

CC_DEPS += \
./doc/examples/main.d 


# Each subdirectory must supply rules for building sources it contributes
doc/examples/%.o: ../doc/examples/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: Cygwin C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


