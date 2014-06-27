/*
    LoopTK: Protein Loop Kinematic Toolkit
    Copyright (C) 2007 Stanford University

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#include "CommandLine.h"
#include <assert.h>
#include <iostream>
using namespace CommandLine;
using namespace std;

namespace CommandLine {

bool string_is_bool(const char* arg, bool& val)
{
	if(arg[0] == '0') val = false;
	else if(arg[0] == '1') val = true;
	else if(0 == strcmp(arg,"true")) val = true;
	else if(0 == strcmp(arg,"false")) val = false;
	else return false;
	return true;
}

bool string_is_int(const char* arg, int& val)
{
	const char* c = arg;
	while(*c) { if(!(isdigit(*c) || (*c)=='-')) return false; c++; }
	val = atoi(arg);
	return true;
}

bool string_is_float(const char* arg, float& val)
{
	const char* c = arg;
	while(*c) { if(!(isdigit(*c) || (*c)=='.' || (*c)=='-')) return false; c++; }
	val = (float)atof(arg);
	return true;
}

bool string_is_double(const char* arg, double& val)
{
	const char* c = arg;
	while(*c) { if(!(isdigit(*c) || (*c)=='.' || (*c)=='-')) return false; c++; }
	val = atof(arg);
	return true;
}


void CommandBase::PrintDescription(std::ostream& out)
{
	out<<"-"<<Name();
	for(size_t i=0; i<args.size(); i++) {
		out << " " << args[i].name;
	}
	out<<": "<<Description()<<endl;
}

void CommandBase::PrintDetailedDescription(ostream& out)
{
	out<<"-"<<Name()<<":"<<endl;
	out<<"Args:"<<endl;
	for(size_t i=0;i<args.size();i++) {
		switch(args[i].type) {
		case Bool:
			out<<"\tbool "<<args[i].name<<" (default "<<(*(bool*)(args[i].data) ? "true":"false")<<")"<<endl;
			break;
		case Int:
			out<<"\tint "<<args[i].name<<" (default "<<(*(int*)(args[i].data))<<")"<<endl;
			break;
		case Float:
			out<<"\tfloat "<<args[i].name<<" (default "<<(*(float*)(args[i].data))<<")"<<endl;
			break;
		case Double:
			out<<"\tdouble "<<args[i].name<<" (default "<<(*(double*)(args[i].data))<<")"<<endl;
			break;
		case String:
			out<<"\tstring "<<args[i].name<<" (default "<<(*(const char**)(args[i].data))<<")"<<endl;
			break;
		default:
			out<<"Error!"<<endl;
			break;
		}
	}
	if(MaxInputs() < 0)
		out<<"Inputs: at least "<<MinInputs()<<endl;
	else if(MaxInputs() == MinInputs())
		out<<"Inputs: "<<MinInputs()<<endl;
	else
		out<<"Inputs: between "<<MinInputs()<<" and "<<MaxInputs()<<endl;
	out<<"Outputs: "<<NumOutputs()<<endl;
	out<<"Description: "; PrintDescription(out);
}
void CommandBase::AddBoolArg(const char* name, bool& val) {
	args.push_back(CommandArg(Bool,name,&val)); }
void CommandBase::AddIntArg(const char* name, int& val) {
	args.push_back(CommandArg(Int,name,&val)); }
void CommandBase::AddFloatArg(const char* name, float& val) {
	args.push_back(CommandArg(Float,name,&val)); }
void CommandBase::AddDoubleArg(const char* name, double& val) {
	args.push_back(CommandArg(Double,name,&val)); }
void CommandBase::AddStringArg(const char* name, const char*& val) {
	args.push_back(CommandArg(String,name,&val)); }

bool CommandBase::ProcessLine()
{
	assert(line != NULL);
	if(line->size() < NumArgs()) {
		cout<<"Not enough arguments to "<<Name()<<endl;
		return false;
	}
	int args_remaining = (int)line->size();
	int j;
	for(j=0; j<NumArgs(); j++,args_remaining--) {
		if(!SetArg(j,(*line)[j])) {
			cout<<"Error reading argument "<<GetArgName(j)<<"="<<(*line)[j]<<" to "<<Name()<<endl;
			return false;
		}
	}
	int num_inputs = args_remaining - NumOutputs();
	if(MaxInputs() >= 0) {
		if(num_inputs > MaxInputs()) {
			cout << "The number of max inputs to -" <<Name()<<" was exceeded."<<endl;
			return false;
		}
	}
	if(num_inputs < MinInputs()) {
		cout << "Not enough inputs given to -" <<Name()<<endl;
		return false;
	}
	inputs.resize(num_inputs);
	for(int i=0;i<num_inputs;i++,j++,args_remaining--) {
		inputs[i]=(*line)[j];
	}
	assert(args_remaining == NumOutputs());
	outputs.resize(NumOutputs());
	for(int i=0;i<NumOutputs();i++,j++,args_remaining--) {
		outputs[i]=(*line)[j];
	}
	assert(j == line->size());
	assert(args_remaining == 0);
	return true;
}

bool CommandBase::SetArg(int arg, const char* val)
{
	if(0 == strcmp(val,"#")) return true;  //default arguments
	switch(args(arg).type) {
	case Bool: return string_is_bool(val,*(bool*)args(arg).data);
	case Int: return string_is_int(val,*(int*)args(arg).data);
	case Float: return string_is_float(val,*(float*)args(arg).data);
	case Double: return string_is_double(val,*(double*)args(arg).data);
	case String: *(const char**)args(arg).data = val;
		return true;
	}
	return false;
}

void PrintUsage(const char* appName)
{
	cout << "USAGE: " << appName << " [COMMAND] args inputs outputs ..." << endl;
	cout << "# as an argument denotes using the default argument." << endl;
}

void PrintCommands(CommandBase* commands [], int num_commands)
{
	cout << "COMMANDS:"<<endl;
	for(int i=0; i<num_commands; i++)
		commands[i]->Print_Description(cout);
	cout<<"For more information, type -help commandname"<<endl;
}

struct CommandHelp : public CommandBase
{
	CommandHelp(CommandBase** commands_input, int num_commands_input)
		:commands(commands_input), num_commands(num_commands_input)
	{
		AddStringArg("command",item_name);
	}
	virtual const char* Name() { return "help"; }
	virtual const char* Description() { return "returns help on a command."; }
	virtual int Do() {
		bool found_command = false;
		for(int k=0; k<num_commands && !found_command; k++) {
			if(0 == strcmp(commands[k]->Name(),item_name)) {
				found_command = true;
				commands[k]->Print_Detailed_Description(cout);
			}
		}
		if(!found_command) {
			cout<<"Unknown command given to help :"<<item_name<<endl;
			return -1;
		}
		return 0;
	}

	CommandBase** commands;
	int num_commands;
	const char* item_name;
};

int DoCommand(CommandBase* c, std::vector<const char*>& inputs)
{
	c->line = &inputs;
	if(!c->ProcessLine()) return -1;
	c->line = NULL;
	int res = c->Do();
	if(res <= 0) return res;
	inputs.erase();
	return 1;
}

int Reader(CommandBase* commands [], int num_commands, int argc, const char** argv)
{
	std::vector<const char*> inputs;
	CommandBase* last_command = NULL;
	CommandHelp help_command(commands,num_commands);

	int command_count = 0;
	bool command_chosen = false;
	for(int i=1; i<argc; i++) {
        double tmp;
		if(argv[i][0] == '-' && !string_is_double(argv[i],tmp)) { //command
            if(string_is_double(argv[i],tmp))
			command_chosen = false;
			//check for help queries first
			if(0 == strcmp("help",&argv[i][1])) {
				if(last_command != NULL) {
					int res = DoCommand(last_command,inputs);
					if(res <= 0) return res;
				}
				last_command = &help_command;
				command_count++;
				command_chosen = true;
			}
			//search through rest of commands
			for(int k=0; k<num_commands && !command_chosen; k++) {
				if(0 == strcmp(commands[k]->Name(),&argv[i][1])) {
					//we've found a command
					if(last_command != NULL) {
						int res = Do_Command(last_command,inputs);
						if(res <= 0) return res;
					}
					last_command = commands[k];
					command_chosen = true;
					command_count++;
				}
			}
			if(!command_chosen) {
				printf("Unknown command %s\n", argv[i]);
				PrintCommands(commands, num_commands);
				return -1;
			}
		}
		else {  //add an argument
			inputs.push_back(argv[i]);
			command_chosen = false;
		}
	}
	if(command_count == 0) {
		if(inputs.size() > 0)
			cout<<"Some stray elements were given on the command line"<<endl;
		PrintCommands(commands, num_commands);
	}
	else {
		assert(last_command != NULL);
		int res = DoCommand(last_command,inputs);
		if(res <= 0) return res;
	}
	return 0;
}

};
