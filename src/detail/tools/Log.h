// Copyright (c) 2012-2013, IGN France.
// Copyright (c) 2012-2024, Oslandia.
// Copyright (c) 2024-2026, SFCGAL team.
// SPDX-License-Identifier: LGPL-2.0-or-later

#ifndef SFCGAL_LOG_H_
#define SFCGAL_LOG_H_

#include "SFCGAL/config.h"

#include <format>
#include <ostream>
#include <source_location>
#include <string>

namespace SFCGAL {

/**
 * [Singleton]Logger class
 *
 * @warning saved_lines and co removed (dangerous for memory and could be done
 * in a LogWriter).
 */
class SFCGAL_API Logger {
public:
  /**
   * destructor
   */
  ~Logger();

  /**
   * log level
   */
  using Level = enum { Debug, Info, Warning, Error, Critical };

  /**
   * singleton accessor
   * @return Pointer to the global Logger instance
   */
  static Logger *
  get();

  /**
   * log a message
   * @param level the log level
   * @param message the message to log
   * @param filename the filename (optional)
   * @param lineNumber the line number in the file (optional)
   */
  void
  log(const Level &level, const std::string &message,
      const std::string &filename = "", const int &lineNumber = -1);

  /**
   * get the current log level
   * @return The current logging level
   */
  const Level &
  logLevel() const;
  /**
   * set the log level
   * @param logLevel The new log level to set
   */
  void
  setLogLevel(const Level &logLevel);

private:
  /**
   * current log level
   */
  Level _logLevel = Warning;
  /**
   * display file position?
   */
  bool _displayFilePosition = true;

  /**
   * private constructor
   */
  Logger(std::ostream &);
  /**
   * no copy constructor
   */
  Logger(const Logger &other);

  std::ostream _out;
};

/**
 * @brief get the logger
 * @return Reference to the global logger instance
 */
SFCGAL_API auto
logger() -> Logger &;

} // namespace SFCGAL

inline void
SFCGAL_LOG(const std::string &level, const std::string &msg,
           std::source_location loc = std::source_location::current())
{
  SFCGAL::Logger::Level lvl = SFCGAL::Logger::Info;
  if (level == "DEBUG") {
    lvl = SFCGAL::Logger::Debug;
  }
  if (level == "WARNING") {
    lvl = SFCGAL::Logger::Warning;
  }
  if (level == "ERROR") {
    lvl = SFCGAL::Logger::Error;
  }
  if (level == "CRITICAL") {
    lvl = SFCGAL::Logger::Critical;
  }

  SFCGAL::Logger::get()->log(lvl, msg, loc.file_name(),
                             static_cast<int>(loc.line()));
}

inline void
LOG_NOTICE(const std::string   &msg,
           std::source_location loc = std::source_location::current())
{
  SFCGAL_LOG("NOTICE", msg, loc);
}

inline void
LOG_ERROR(const std::string   &msg,
          std::source_location loc = std::source_location::current())
{
  SFCGAL_LOG("ERROR", msg, loc);
}

#ifndef NDEBUG
inline void
LOG_DEBUG(const std::string   &msg,
          std::source_location loc = std::source_location::current())
{
  SFCGAL_LOG("DEBUG", msg, loc);
}
#else
inline void
LOG_DEBUG(const std::string & /*msg*/,
          std::source_location /*loc*/ = std::source_location::current())
{
  // no-op in release
}
#endif

/**
 *
 * Helper method to log debug message
 *
 * \code
 * SFCGAL_DEBUG( "start new method" ) ;
 * \endcode
 */
inline void
SFCGAL_DEBUG(const std::string   &message,
             std::source_location loc = std::source_location::current())
{
  SFCGAL::Logger::get()->log(SFCGAL::Logger::Debug, message, loc.file_name(),
                             static_cast<int>(loc.line()));
}

/**
 *
 * Helper method to log information message
 *
 * \code
 * SFCGAL_INFO( "start new method" ) ;
 * \endcode
 */
inline void
SFCGAL_INFO(const std::string   &message,
            std::source_location loc = std::source_location::current())
{
  SFCGAL::Logger::get()->log(SFCGAL::Logger::Info, message, loc.file_name(),
                             static_cast<int>(loc.line()));
}

/**
 *
 * Helper method to log warning message
 *
 * \code
 * SFCGAL_WARNING( "start new method" ) ;
 * \endcode
 */
inline void
SFCGAL_WARNING(const std::string   &message,
               std::source_location loc = std::source_location::current())
{
  SFCGAL::Logger::get()->log(SFCGAL::Logger::Warning, message, loc.file_name(),
                             static_cast<int>(loc.line()));
}

/**
 *
 * Helper method to log error message
 *
 * \code
 * SFCGAL_ERROR( "invalid geometry" ) ;
 * \endcode
 */
inline void
SFCGAL_ERROR(const std::string   &message,
             std::source_location loc = std::source_location::current())
{
  SFCGAL::Logger::get()->log(SFCGAL::Logger::Error, message, loc.file_name(),
                             static_cast<int>(loc.line()));
}

/**
 *
 * Helper method to log critical message
 *
 * \code
 * SFCGAL_ERROR( "unexpected behavior in triangulate" ) ;
 * \endcode
 */
inline void
SFCGAL_CRITICAL(const std::string   &message,
                std::source_location loc = std::source_location::current())
{
  SFCGAL::Logger::get()->log(SFCGAL::Logger::Critical, message, loc.file_name(),
                             static_cast<int>(loc.line()));
}

#endif // SFCGAL_LOG_H_
