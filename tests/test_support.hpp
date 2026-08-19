#pragma once

#include <sstream>
#include <stdexcept>
#include <string>

namespace wcns::test {

inline void require(bool condition, const char* expression, const char* file, int line)
{
    if (condition) {
        return;
    }
    std::ostringstream message;
    message << file << ':' << line << ": requirement failed: " << expression;
    throw std::runtime_error(message.str());
}

template<class Exception, class Function>
void require_throws(Function&& function, const char* expression, const char* file, int line)
{
    try {
        function();
    } catch (const Exception&) {
        return;
    } catch (...) {
        std::ostringstream message;
        message << file << ':' << line << ": wrong exception from: " << expression;
        throw std::runtime_error(message.str());
    }
    std::ostringstream message;
    message << file << ':' << line << ": expected exception from: " << expression;
    throw std::runtime_error(message.str());
}

} // namespace wcns::test

#define WCNS_REQUIRE(expression) \
    ::wcns::test::require(static_cast<bool>(expression), #expression, __FILE__, __LINE__)

#define WCNS_REQUIRE_THROWS(exception_type, expression) \
    ::wcns::test::require_throws<exception_type>( \
        [&] { static_cast<void>(expression); }, #expression, __FILE__, __LINE__)

