#ifndef tt_message_h
#define tt_message_h

#include <string>

class TtMessage;

namespace tt
{
	extern TtMessage* g_msg;
	extern int g_verbose_level;
}

enum TtVerboseLevel
{
	TT_SILENT_MODE,
	TT_BATCH_MODE,
	TT_NORMAL_MODE,
	TT_VERBOSE_MODE,
};

class TtMessage
{
public:
	TtMessage();
	virtual ~TtMessage();

	virtual void printError(const char* str);
	virtual void printInfo(const char* str);
	virtual void printVerbose(const char* str);
	virtual void printData(const char* str);

	virtual void startProgress(int percent);
	virtual void printProgress(int percent);
	virtual void endProgress(int percent);

	std::string m_history;
};

int  tt_get_verbose_level();
void tt_set_verbose_level(int lv);
void tt_set_message_class(TtMessage* o);

std::string tt_get_history();
void tt_error(const char* fmt, ...);
void tt_info(const char* fmt, ...);
void tt_verbose(const char* fmt, ...);
void tt_print_data(const char* fmt, ...);
void tt_print_text(const char** text);

void tt_progress_start(int percent);
void tt_progress(int percent);
void tt_progress_end(int percent);

#endif  // tt_message_h

