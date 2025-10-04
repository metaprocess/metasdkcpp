#ifndef COMMANDEXECUTOR_H
#define COMMANDEXECUTOR_H

#include <string>

class CommandExecutor
{
public:
    CommandExecutor();

	static std::string runCommand(const std::string& cmd);
};

#endif // COMMANDEXECUTOR_H
