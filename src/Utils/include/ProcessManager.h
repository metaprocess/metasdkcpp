#ifndef PROCESSMANAGER_H
#define PROCESSMANAGER_H

#include <iostream>
#include <vector>
#include <string>
#include <array>
#include <unistd.h>
#include <sys/wait.h>
#include <fcntl.h>
#include <thread>
#include <functional>
#include <map>

class ProcessManager {
public:
    ProcessManager();
    ~ProcessManager();

    pid_t startProcess(const std::string& command, const std::vector<std::string>& args, const std::map<std::string, std::string>& env, std::function<void(const std::string&)> outputCallback);
    void writeToProcess(const std::string& data);
    void killProcess();

private:
    void readFromProcess(std::function<void(const std::string&)> outputCallback);

    pid_t pid;
    int stdin_fd;
    int stdout_fd;
    int stderr_fd;
    // std::thread read_thread;
};

#endif // PROCESSMANAGER_H