#include <stdio.h>
#include <stdarg.h>
#include "message.h"
#include "string.h"

#define BUFFER_SIZE 1024

using namespace std;

namespace tt
{
	TtMessage* g_msg = 0;
	int g_verbose_level = TT_NORMAL_MODE;
}

class InitMessage
{
public:
	InitMessage()
	{
		tt::g_msg = new TtMessage();
	}
};

InitMessage init_msg;

//=========================================================================

TtMessage::TtMessage()
{
}

TtMessage::~TtMessage()
{
}

void
TtMessage::printError(const char* str)
{
	fprintf(stderr, "%s", str);
	fflush(stderr);
}

void
TtMessage::printInfo(const char* str)
{
	fprintf(stdout, "%s", str);
	fflush(stdout);

	m_history = str;
	tt_remove_newline(m_history);
}

void
TtMessage::printVerbose(const char* str)
{
	fprintf(stdout, "%s", str);
	fflush(stdout);

	m_history = str;
	tt_remove_newline(m_history);
}

void
TtMessage::printData(const char* str)
{
	fprintf(stdout, "%s", str);
	fflush(stdout);
}

void
TtMessage::startProgress(int percent)
{
	fprintf(stderr, "[%3d %%]", percent);
	fflush(stderr);
}

void
TtMessage::printProgress(int percent)
{
	fprintf(stderr, "\b\b\b\b\b\b\b[%3d %%]", percent);
	fflush(stderr);
}

void
TtMessage::endProgress(int percent)
{
	fprintf(stderr, "\b\b\b\b\b\b\b");
	fflush(stderr);
}

//=========================================================================

int
tt_get_verbose_level()
{
	return tt::g_verbose_level;
}

void
tt_set_verbose_level(int lv)
{
	tt::g_verbose_level = lv;
}

void
tt_set_message_class(TtMessage* o)
{
	if (tt::g_msg)
	{
		delete tt::g_msg;
		tt::g_msg = 0;
	}
	tt::g_msg = o;
}

//=========================================================================

std::string
tt_get_history()
{
	return tt::g_msg->m_history;
}

void
tt_error(const char* fmt, ...)
{
	if (tt::g_verbose_level < TT_BATCH_MODE) return;

	char buffer[BUFFER_SIZE];
	va_list ap;

	va_start(ap, fmt);
	vsprintf(buffer, fmt, ap);
	va_end(ap);

	tt::g_msg->printError(buffer);
}

void
tt_info(const char* fmt, ...)
{
	if (tt::g_verbose_level < TT_BATCH_MODE) return;

	char buffer[BUFFER_SIZE];
	va_list ap;

	va_start(ap, fmt);
	vsprintf(buffer, fmt, ap);
	va_end(ap);

	tt::g_msg->printInfo(buffer);
}

void
tt_verbose(const char* fmt, ...)
{
	if (tt::g_verbose_level < TT_VERBOSE_MODE) return;

	char buffer[BUFFER_SIZE];
	va_list ap;

	va_start(ap, fmt);
	vsprintf(buffer, fmt, ap);
	va_end(ap);

	tt::g_msg->printVerbose(buffer);
}

void
tt_print_data(const char* fmt, ...)
{
	if (tt::g_verbose_level < TT_BATCH_MODE) return;

	char buffer[BUFFER_SIZE];
	va_list ap;

	va_start(ap, fmt);
	vsprintf(buffer, fmt, ap);
	va_end(ap);

	tt::g_msg->printData(buffer);
}

void
tt_print_text(const char** text)
{
	if (tt::g_verbose_level < TT_BATCH_MODE) return;

	for (int i=0; text[i][0]; i++)
	{
		tt::g_msg->printInfo(text[i]);
	}
}

void
tt_progress_start(int percent)
{
	if (tt::g_verbose_level < TT_BATCH_MODE) return;

	tt::g_msg->startProgress(percent);
}

void
tt_progress(int percent)
{
	if (tt::g_verbose_level < TT_BATCH_MODE) return;

	tt::g_msg->printProgress(percent);
}

void
tt_progress_end(int percent)
{
	if (tt::g_verbose_level < TT_BATCH_MODE) return;

	tt::g_msg->endProgress(percent);
}

