#include "tools/jasonParser.h"

#include <sstream>
#include <fstream>


using namespace std;
using namespace jParser;

string deserialize(const string& ref) {
    string out = "";
    for (size_t i = 0; i < ref.length(); i++) {
        if (ref[i] == '\\' && i + 1 < ref.length()) {
            int plus = 2;
            if (ref[i + 1] == '\"') {
                out += '"';
            }
            else if (ref[i + 1] == '\\') {
                out += '\\';
            }
            else if (ref[i + 1] == '/') {
                out += '/';
            }
            else if (ref[i + 1] == 'b') {
                out += '\b';
            }
            else if (ref[i + 1] == 'f') {
                out += '\f';
            }
            else if (ref[i + 1] == 'n') {
                out += '\n';
            }
            else if (ref[i + 1] == 'r') {
                out += '\r';
            }
            else if (ref[i + 1] == 't') {
                out += '\t';
            }
            else if (ref[i + 1] == 'u' && i + 5 < ref.length()) {
                unsigned long long v = 0;
                for (int j = 0; j < 4; j++) {
                    v *= 16;
                    if (ref[i + 2 + j] <= '9' && ref[i + 2 + j] >= '0') v += ref[i + 2 + j] - '0';
                    if (ref[i + 2 + j] <= 'f' && ref[i + 2 + j] >= 'a') v += ref[i + 2 + j] - 'a' + 10;
                }
                out += (char)v;
                plus = 6;
            }
            i += plus - 1;
            continue;
        }
        out += ref[i];
    }
    return out;
}

string jValue::makesp(int d) const {
    string s = "";
    while (d--) s += "  ";
    return s;
}
string jValue::to_string_d(int d) const {
    if (type == jType::JSTRING)   return string("\"") + svalue + string("\"");
    if (type == jType::JNUMBER)   return svalue;
    if (type == jType::JBOOLEAN)  return svalue;
    if (type == jType::JNULL)     return "null";
    if (type == jType::JOBJECT) {
        string s = string("{\n");
        for (size_t i = 0; i < properties.size(); i++) {
            s += makesp(d) + string("\"") + properties[i].first + string("\": ") + properties[i].second.to_string_d(d + 1) + string(i == properties.size() - 1 ? "" : ",") + string("\n");
        }
        s += makesp(d - 1) + string("}");
        return s;
    }
    if (type == jType::JARRAY) {
        string s = "[";
        for (size_t i = 0; i < arr.size(); i++) {
            if (i) s += ", ";
            s += arr[i].to_string_d(d + 1);
        }
        s += "]";
        return s;
    }
    return "##";
}
jValue::jValue() {
    this->type = jType::JUNKNOWN;
}
jValue::jValue(jType tp) {
    this->type = tp;
}

string jValue::to_string() const {
    return to_string_d(1);
}
jType jValue::get_type() const {
    return type;
}
void jValue::set_type(jType tp) {
    type = tp;
}
void jValue::add_property(string key, jValue v) {
    mpindex[key] = properties.size();
    properties.push_back(make_pair(key, v));
}
void jValue::add_element(jValue v) {
    arr.push_back(v);
}
void jValue::set_string(string s) {
    svalue = s;
}
int jValue::as_int() const {
    stringstream ss;
    ss << svalue;
    int k;
    ss >> k;
    return k;
}
float jValue::as_float() const {
    stringstream ss;
    ss << svalue;
    float k;
    ss >> k;
    return k;
}
double jValue::as_double() const {
    stringstream ss;
    ss << svalue;
    double k;
    ss >> k;
    return k;
}
bool jValue::as_bool() const {
    if (svalue == "true") return true;
    return false;
}
void* jValue::as_null() const {
    return NULL;
}
string jValue::as_string() const {
    return deserialize(svalue);
}
int jValue::size() const {
    if (type == jType::JARRAY) {
        return (int)arr.size();
    }
    if (type == jType::JOBJECT) {
        return (int)properties.size();;
    }
    return 0;
}
bool jValue::contains(const std::string& s) const {
    return mpindex.find(s) != mpindex.end();
}
std::string jValue::getKey(int i) const {
    if (type == jType::JARRAY) {
        return std::to_string(i);
    }
    if (type == jType::JOBJECT) {
        return properties[i].first;
    }
}
jValue jValue::operator[](int i) {
    if (type == jType::JARRAY) {
        return arr[i];
    }
    if (type == jType::JOBJECT) {
        return properties[i].second;
    }
    return jValue();
}
jValue jValue::operator[](std::string s) {
    if (mpindex.find(s) == mpindex.end()) return jValue();
    return properties[mpindex[s]].second;
}

struct parser::token {
    string value;
    token_type type;
    token(string value = "", token_type type = token_type::UNKNOWN) : value(value), type(type) {}
};
bool parser::is_whitespace(const char c) {
    return isspace(c);
}
int parser::next_whitespace(const string& source, int i) {
    while (i < (int)source.length()) {
        if (source[i] == '"') {
            i++;
            while (i < (int)source.length() && (source[i] != '"' || source[i - 1] == '\\')) i++;
        }
        if (source[i] == '\'') {
            i++;
            while (i < (int)source.length() && (source[i] != '\'' || source[i - 1] == '\\')) i++;
        }
        if (is_whitespace(source[i])) return i;
        i++;
    }
    return (int)source.length();
}
int parser::skip_whitespaces(const string& source, int i) {
    while (i < (int)source.length()) {
        if (!is_whitespace(source[i])) return i;
        i++;
    }
    return -1;
}

vector<parser::token> parser::tokenize(string source) {
    source += " ";
    vector<token> tokens;
    int index = skip_whitespaces(source, 0);
    while (index >= 0) {
        int next = next_whitespace(source, index);
        string str = source.substr(index, next - index);

        size_t k = 0;
        while (k < str.length()) {
            if (str[k] == '"') {
                size_t tmp_k = k + 1;
                while (tmp_k < str.length() && (str[tmp_k] != '"' || str[tmp_k - 1] == '\\')) tmp_k++;
                tokens.push_back(token(str.substr(k + 1, tmp_k - k - 1), token_type::STRING));
                k = tmp_k + 1;
                continue;
            }
            if (str[k] == '\'') {
                size_t tmp_k = k + 1;
                while (tmp_k < str.length() && (str[tmp_k] != '\'' || str[tmp_k - 1] == '\\')) tmp_k++;
                tokens.push_back(token(str.substr(k + 1, tmp_k - k - 1), token_type::STRING));
                k = tmp_k + 1;
                continue;
            }
            if (str[k] == ',') {
                tokens.push_back(token(",", token_type::COMMA));
                k++;
                continue;
            }
            if (str[k] == 't' && k + 3 < str.length() && str.substr(k, 4) == "true") {
                tokens.push_back(token("true", token_type::BOOLEAN));
                k += 4;
                continue;
            }
            if (str[k] == 'f' && k + 4 < str.length() && str.substr(k, 5) == "false") {
                tokens.push_back(token("false", token_type::BOOLEAN));
                k += 5;
                continue;
            }
            if (str[k] == 'n' && k + 3 < str.length() && str.substr(k, 4) == "null") {
                tokens.push_back(token("null", token_type::NUL));
                k += 4;
                continue;
            }
            if (str[k] == '}') {
                tokens.push_back(token("}", token_type::CROUSH_CLOSE));
                k++;
                continue;
            }
            if (str[k] == '{') {
                tokens.push_back(token("{", token_type::CROUSH_OPEN));
                k++;
                continue;
            }
            if (str[k] == ']') {
                tokens.push_back(token("]", token_type::BRACKET_CLOSE));
                k++;
                continue;
            }
            if (str[k] == '[') {
                tokens.push_back(token("[", token_type::BRACKET_OPEN));
                k++;
                continue;
            }
            if (str[k] == ':') {
                tokens.push_back(token(":", token_type::COLON));
                k++;
                continue;
            }
            if (str[k] == '-' || (str[k] <= '9' && str[k] >= '0')) {
                size_t tmp_k = k;
                if (str[tmp_k] == '-') tmp_k++;
                while (tmp_k < str.size() && ((str[tmp_k] <= '9' && str[tmp_k] >= '0') || str[tmp_k] == '.')) tmp_k++;
                tokens.push_back(token(str.substr(k, tmp_k - k), token_type::NUMBER));
                k = tmp_k;
                continue;
            }
            tokens.push_back(token(str.substr(k), token_type::UNKNOWN));
            k = str.length();
        }

        index = skip_whitespaces(source, next);
    }
    // for (int i=0;i<tokens.size();i++) {
      // cout << i << " " << tokens[i].value << endl;
    // }
    return tokens;
}


jValue parser::json_parse(vector<token> v, int i, int& r) {
    jValue current;
    if (v[i].type == token_type::CROUSH_OPEN) {
        current.set_type(jType::JOBJECT);
        int k = i + 1;
        while (v[k].type != token_type::CROUSH_CLOSE) {
            string key = v[k].value;
            k += 2; // k+1 should be ':'
            int j = k;
            jValue vv = json_parse(v, k, j);
            current.add_property(key, vv);
            k = j;
            if (v[k].type == token_type::COMMA) k++;
        }
        r = k + 1;
        return current;
    }
    if (v[i].type == token_type::BRACKET_OPEN) {
        current.set_type(jType::JARRAY);
        int k = i + 1;
        while (v[k].type != token_type::BRACKET_CLOSE) {
            int j = k;
            jValue vv = json_parse(v, k, j);
            current.add_element(vv);
            k = j;
            if (v[k].type == token_type::COMMA) k++;
        }
        r = k + 1;
        return current;
    }
    if (v[i].type == token_type::NUMBER) {
        current.set_type(jType::JNUMBER);
        current.set_string(v[i].value);
        r = i + 1;
        return current;
    }
    if (v[i].type == token_type::STRING) {
        current.set_type(jType::JSTRING);
        current.set_string(v[i].value);
        r = i + 1;
        return current;
    }
    if (v[i].type == token_type::BOOLEAN) {
        current.set_type(jType::JBOOLEAN);
        current.set_string(v[i].value);
        r = i + 1;
        return current;
    }
    if (v[i].type == token_type::NUL) {
        current.set_type(jType::JNULL);
        current.set_string("null");
        r = i + 1;
        return current;
    }
    return current;
}

jValue parser::parse(const string& str) {
    int k;
    return json_parse(tokenize(str), 0, k);
}
jValue parser::parse_file(const string& filename) {
    ifstream in(filename.c_str());
    string str = "";
    string tmp;
    while (getline(in, tmp)) str += tmp;
    in.close();
    return parser::parse(str);
}