#include "ProcessManager.h"
#include <future>
#include <Utils.h>

ProcessManager::ProcessManager() : pid(-1), stdin_fd(-1), stdout_fd(-1), stderr_fd(-1) {}

ProcessManager::~ProcessManager() {
    killProcess();
}

pid_t ProcessManager::startProcess(const std::string& command, const std::vector<std::string>& args, const std::map<std::string, std::string>& env, std::function<void(const std::string&)> outputCallback) {
    pid_t child_pid;
    int pipe_stdin[2], pipe_stdout[2], pipe_stderr[2];

    if (pipe(pipe_stdin) == -1 || pipe(pipe_stdout) == -1 || pipe(pipe_stderr) == -1) {
        perror("pipe");
        return -1;
    }

    child_pid = fork();

    if (child_pid == -1) {
        perror("fork");
        return -1;
    } else if (child_pid == 0) {
        // Child process
        close(pipe_stdin[1]);  // Close write end of stdin pipe
        close(pipe_stdout[0]); // Close read end of stdout pipe
        close(pipe_stderr[0]); // Close read end of stderr pipe

        dup2(pipe_stdin[0], STDIN_FILENO);  // Redirect stdin
        dup2(pipe_stdout[1], STDOUT_FILENO); // Redirect stdout
        dup2(pipe_stderr[1], STDERR_FILENO); // Redirect stdout

        // Set environment variables
        for (const auto& pair : env) {
            setenv(pair.first.c_str(), pair.second.c_str(), 1);
        }

        std::vector<char*> c_args;
        c_args.push_back(const_cast<char*>(command.c_str()));
        for (const auto& arg : args) {
            c_args.push_back(const_cast<char*>(arg.c_str()));
        }
        c_args.push_back(nullptr);

        execvp(command.c_str(), c_args.data());
        perror("execvp");
        _exit(1);
    } else {
        // Parent process
        close(pipe_stdin[0]);  // Close read end of stdin pipe
        close(pipe_stdout[1]); // Close write end of stdout pipe
        close(pipe_stderr[1]); // Close write end of stderr pipe

        stdin_fd = pipe_stdin[1];
        stdout_fd = pipe_stdout[0];
        stderr_fd = pipe_stderr[0];
    }

    pid = child_pid;

    // Start a thread to read from stdout and call the callback
    std::thread* read_thread{nullptr};
    read_thread = new std::thread([&command, this, outputCallback, &read_thread]() {
        auto _str = command.substr(0, 16);
        Utils::setThreadName(_str.c_str());
        readFromProcess(outputCallback);
        killProcess();
        // delete read_thread;
        read_thread = nullptr;
        // std::cerr << "exiting Process manager read loop\n";
    });

    return child_pid;
}
void ProcessManager::writeToProcess(const std::string& data) {
    if (stdin_fd == -1) {
        std::cerr << "Process not started or stdin not available" << std::endl;
        return;
    }
    write(stdin_fd, data.c_str(), data.size());
}

void ProcessManager::killProcess() {
    if (pid > 0) {
        kill(pid, SIGTERM);
        waitpid(pid, nullptr, 0);
        close(stdin_fd);
        close(stdout_fd);
        close(stderr_fd);
        pid = -1;
        stdin_fd = -1;
        stdout_fd = -1;
        stderr_fd = -1;

        // if (read_thread.joinable()) {
        //     read_thread.join();
        // }
    }
}

void ProcessManager::readFromProcess(std::function<void(const std::string&)> outputCallback) {
    if (stdout_fd == -1 || stderr_fd == -1) {
        std::cerr << "Process not started or stdout/stderr not available" << std::endl;
        return;
    }

    std::array<char, 128> buffer;
    ssize_t bytes_read;
    fd_set read_fds;
    int max_fd = std::max(stdout_fd, stderr_fd);

    while (true) {
        FD_ZERO(&read_fds);
        FD_SET(stdout_fd, &read_fds);
        FD_SET(stderr_fd, &read_fds);

        if (select(max_fd + 1, &read_fds, nullptr, nullptr, nullptr) == -1) {
            perror("select");
            break;
        }

        if (FD_ISSET(stdout_fd, &read_fds)) {
            bytes_read = read(stdout_fd, buffer.data(), buffer.size());
            if (bytes_read > 0) {
                std::string result(buffer.data(), bytes_read);
                if(outputCallback)
                {
                    outputCallback(result);
                }
            } else if (bytes_read == 0) {
                break; // EOF
            }
        }

        if (FD_ISSET(stderr_fd, &read_fds)) {
            bytes_read = read(stderr_fd, buffer.data(), buffer.size());
            if (bytes_read > 0) {
                std::string result(buffer.data(), bytes_read);
                if(outputCallback)
                {
                    outputCallback(result);
                }
            } else if (bytes_read == 0) {
                break; // EOF
            }
        }
    }

}
