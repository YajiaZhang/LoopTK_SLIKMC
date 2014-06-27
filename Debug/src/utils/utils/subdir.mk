################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/utils/utils/CommandLine.cpp \
../src/utils/utils/ProgressPrinter.cpp \
../src/utils/utils/RefPointer.cpp \
../src/utils/utils/SignalHandler.cpp \
../src/utils/utils/SimpleParser.cpp \
../src/utils/utils/Trace.cpp \
../src/utils/utils/fileutils.cpp \
../src/utils/utils/ioutils.cpp \
../src/utils/utils/set.cpp \
../src/utils/utils/stringutils.cpp \
../src/utils/utils/unionfind.cpp 

OBJS += \
./src/utils/utils/CommandLine.o \
./src/utils/utils/ProgressPrinter.o \
./src/utils/utils/RefPointer.o \
./src/utils/utils/SignalHandler.o \
./src/utils/utils/SimpleParser.o \
./src/utils/utils/Trace.o \
./src/utils/utils/fileutils.o \
./src/utils/utils/ioutils.o \
./src/utils/utils/set.o \
./src/utils/utils/stringutils.o \
./src/utils/utils/unionfind.o 

CPP_DEPS += \
./src/utils/utils/CommandLine.d \
./src/utils/utils/ProgressPrinter.d \
./src/utils/utils/RefPointer.d \
./src/utils/utils/SignalHandler.d \
./src/utils/utils/SimpleParser.d \
./src/utils/utils/Trace.d \
./src/utils/utils/fileutils.d \
./src/utils/utils/ioutils.d \
./src/utils/utils/set.d \
./src/utils/utils/stringutils.d \
./src/utils/utils/unionfind.d 


# Each subdirectory must supply rules for building sources it contributes
src/utils/utils/%.o: ../src/utils/utils/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cygwin C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


