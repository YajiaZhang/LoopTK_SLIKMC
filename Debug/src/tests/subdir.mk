################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../src/tests/SeedSampling.cc \
../src/tests/atom.cc \
../src/tests/atomshell.cc \
../src/tests/basic.cc \
../src/tests/blockbondshell.cc \
../src/tests/collision.cc \
../src/tests/destruct.cc \
../src/tests/init.cc \
../src/tests/math.cc \
../src/tests/pdbio.cc \
../src/tests/pdbio2.cc \
../src/tests/pdbio3.cc \
../src/tests/phipsi.cc \
../src/tests/rotation.cc \
../src/tests/utilities.cc 

OBJS += \
./src/tests/SeedSampling.o \
./src/tests/atom.o \
./src/tests/atomshell.o \
./src/tests/basic.o \
./src/tests/blockbondshell.o \
./src/tests/collision.o \
./src/tests/destruct.o \
./src/tests/init.o \
./src/tests/math.o \
./src/tests/pdbio.o \
./src/tests/pdbio2.o \
./src/tests/pdbio3.o \
./src/tests/phipsi.o \
./src/tests/rotation.o \
./src/tests/utilities.o 

CC_DEPS += \
./src/tests/SeedSampling.d \
./src/tests/atom.d \
./src/tests/atomshell.d \
./src/tests/basic.d \
./src/tests/blockbondshell.d \
./src/tests/collision.d \
./src/tests/destruct.d \
./src/tests/init.d \
./src/tests/math.d \
./src/tests/pdbio.d \
./src/tests/pdbio2.d \
./src/tests/pdbio3.d \
./src/tests/phipsi.d \
./src/tests/rotation.d \
./src/tests/utilities.d 


# Each subdirectory must supply rules for building sources it contributes
src/tests/%.o: ../src/tests/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: Cygwin C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


