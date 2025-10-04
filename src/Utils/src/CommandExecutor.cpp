#include "CommandExecutor.h"

#include <stdio.h>
#include "Definitions.h"

CommandExecutor::CommandExecutor()
{
//xrandr | sed -n  's|\(^.*\) connected .*|\1|p'
    //    route -n| grep  '^0\.0\.0\.0' | tr -s ' ' ' ' | cut -d' ' -f8
}

bool endsWith(const std::string& str, const std::string& suffix) {
    if (str.length() < suffix.length()) {
        return false;
    }
    return str.compare(str.length() - suffix.length(), suffix.length(), suffix) == 0;
}

std::string CommandExecutor::runCommand(const std::string &cmd)
{
    std::string outStr;
    FILE* fid = popen(cmd.c_str(), "r");
    const int len = 1e3;
    auto buffer = std::make_unique<char[]>(len);
    if(fid)
    {
        while (!feof(fid))
        {
            if (fgets(buffer.get(), len, fid) != nullptr)
            {
               outStr.append(buffer.get());
            }
        }
        pclose(fid);
    }
    if(endsWith(outStr, "\n"))
    {
        outStr.erase(outStr.length() - 1);
    }
    return std::move(outStr);
}
