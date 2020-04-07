#include <Windows.h>
#include <stdexcept>
#include <string>
#include <vector>

struct Process_Data {
	LPSTR commandline;
	STARTUPINFO startupinfo;
	PROCESS_INFORMATION process_information;
};

void plot_result()
{
	std::vector<Process_Data> list_of_process_data(1);

	list_of_process_data[0].commandline = R"(C:\Program Files\gnuplot\bin\gnuplot.exe plot.gp -persist)";


	for (auto & process_data : list_of_process_data)
	{
		if (!CreateProcess(nullptr,
			process_data.commandline,
			nullptr,
			nullptr,
			false,
			0,
			nullptr,
			nullptr,
			&process_data.startupinfo,
			&process_data.process_information))
		{
			throw std::runtime_error("Failed to create process (Error code " + std::to_string(GetLastError()) + ")");
		}
	}

	
	for (auto & process_data : list_of_process_data)
	{
		WaitForSingleObject(process_data.process_information.hProcess, INFINITE);
		CloseHandle(process_data.process_information.hProcess);
		CloseHandle(process_data.process_information.hThread);
	}
}