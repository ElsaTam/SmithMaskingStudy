/**
https://github.com/amir-s/jute/blob/master/jute.h
*/

#pragma once

#include <map>
#include <string>
#include <vector>

namespace jParser {

    enum class jType { JSTRING, JOBJECT, JARRAY, JBOOLEAN, JNUMBER, JNULL, JUNKNOWN };

    class jValue {
    private:
        std::string svalue;
        jType type;
        std::vector<std::pair<std::string, jValue> > properties;
        std::map<std::string, size_t> mpindex;
        std::vector<jValue> arr;
        std::string makesp(int) const;
        std::string to_string_d(int) const ;
    public:
        jValue();
        jValue(jType);

        std::string to_string() const;
        jType get_type() const;
        void set_type(jType);
        void add_property(std::string key, jValue v);
        void add_element(jValue v);
        void set_string(std::string s);
        int as_int() const;
        float as_float() const;
        double as_double() const;
        bool as_bool() const;
        void* as_null() const;
        std::string as_string() const;
        int size() const;
        bool contains(const std::string& s) const;
        std::string getKey(int i) const;
        jValue operator[](int i);
        jValue operator[](std::string s);
    };

    class parser {
    private:
        enum class token_type { UNKNOWN, STRING, NUMBER, CROUSH_OPEN, CROUSH_CLOSE, BRACKET_OPEN, BRACKET_CLOSE, COMMA, COLON, BOOLEAN, NUL };

        struct token;
        static bool is_whitespace(const char c);
        static int next_whitespace(const std::string& source, int i);
        static int skip_whitespaces(const std::string& source, int i);

        static std::vector<token> tokenize(std::string source);
        static jValue json_parse(std::vector<token> v, int i, int& r);
    public:
        static jValue parse(const std::string& str);
        static jValue parse_file(const std::string& str);
    };
}