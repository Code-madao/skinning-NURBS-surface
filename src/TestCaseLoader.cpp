#include "../include/TestCaseLoader.h"

#include <cctype>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <variant>

#include <glm/glm.hpp>

namespace
{
    struct JsonValue;

    using JsonArray = std::vector<JsonValue>;
    using JsonObject = std::map<std::string, JsonValue>;

    struct JsonValue
    {
        using Storage = std::variant<std::nullptr_t, double, std::string, JsonArray, JsonObject>;
        Storage value;

        bool isNumber() const { return std::holds_alternative<double>(value); }
        bool isString() const { return std::holds_alternative<std::string>(value); }
        bool isArray() const { return std::holds_alternative<JsonArray>(value); }
        bool isObject() const { return std::holds_alternative<JsonObject>(value); }

        double asNumber() const
        {
            if (!isNumber())
                throw std::runtime_error("Expected JSON number.");
            return std::get<double>(value);
        }

        const std::string &asString() const
        {
            if (!isString())
                throw std::runtime_error("Expected JSON string.");
            return std::get<std::string>(value);
        }

        const JsonArray &asArray() const
        {
            if (!isArray())
                throw std::runtime_error("Expected JSON array.");
            return std::get<JsonArray>(value);
        }

        const JsonObject &asObject() const
        {
            if (!isObject())
                throw std::runtime_error("Expected JSON object.");
            return std::get<JsonObject>(value);
        }
    };

    class JsonParser
    {
    public:
        explicit JsonParser(std::string text) : text_(std::move(text)), pos_(0) {}

        JsonValue parse()
        {
            JsonValue value = parseValue();
            skipWhitespace();
            if (pos_ != text_.size())
                throw std::runtime_error("Unexpected trailing characters in JSON.");
            return value;
        }

    private:
        JsonValue parseValue()
        {
            skipWhitespace();
            if (pos_ >= text_.size())
                throw std::runtime_error("Unexpected end of JSON.");

            char ch = text_[pos_];
            if (ch == '{')
                return JsonValue{parseObject()};
            if (ch == '[')
                return JsonValue{parseArray()};
            if (ch == '"')
                return JsonValue{parseString()};
            if (ch == '-' || std::isdigit(static_cast<unsigned char>(ch)))
                return JsonValue{parseNumber()};
            if (matchLiteral("null"))
                return JsonValue{nullptr};

            throw std::runtime_error("Unsupported JSON token.");
        }

        JsonObject parseObject()
        {
            expect('{');
            JsonObject object;
            skipWhitespace();
            if (peek('}'))
            {
                expect('}');
                return object;
            }

            while (true)
            {
                skipWhitespace();
                std::string key = parseString();
                skipWhitespace();
                expect(':');
                object.emplace(std::move(key), parseValue());
                skipWhitespace();
                if (peek('}'))
                {
                    expect('}');
                    break;
                }
                expect(',');
            }
            return object;
        }

        JsonArray parseArray()
        {
            expect('[');
            JsonArray array;
            skipWhitespace();
            if (peek(']'))
            {
                expect(']');
                return array;
            }

            while (true)
            {
                array.push_back(parseValue());
                skipWhitespace();
                if (peek(']'))
                {
                    expect(']');
                    break;
                }
                expect(',');
            }
            return array;
        }

        std::string parseString()
        {
            expect('"');
            std::string result;
            while (pos_ < text_.size())
            {
                char ch = text_[pos_++];
                if (ch == '"')
                    return result;
                if (ch == '\\')
                {
                    if (pos_ >= text_.size())
                        throw std::runtime_error("Invalid JSON escape.");
                    char escaped = text_[pos_++];
                    switch (escaped)
                    {
                    case '"':
                    case '\\':
                    case '/':
                        result.push_back(escaped);
                        break;
                    case 'b':
                        result.push_back('\b');
                        break;
                    case 'f':
                        result.push_back('\f');
                        break;
                    case 'n':
                        result.push_back('\n');
                        break;
                    case 'r':
                        result.push_back('\r');
                        break;
                    case 't':
                        result.push_back('\t');
                        break;
                    default:
                        throw std::runtime_error("Unsupported JSON escape sequence.");
                    }
                    continue;
                }
                result.push_back(ch);
            }
            throw std::runtime_error("Unterminated JSON string.");
        }

        double parseNumber()
        {
            size_t start = pos_;
            if (text_[pos_] == '-')
                ++pos_;
            while (pos_ < text_.size() && std::isdigit(static_cast<unsigned char>(text_[pos_])))
                ++pos_;
            if (pos_ < text_.size() && text_[pos_] == '.')
            {
                ++pos_;
                while (pos_ < text_.size() && std::isdigit(static_cast<unsigned char>(text_[pos_])))
                    ++pos_;
            }
            if (pos_ < text_.size() && (text_[pos_] == 'e' || text_[pos_] == 'E'))
            {
                ++pos_;
                if (pos_ < text_.size() && (text_[pos_] == '+' || text_[pos_] == '-'))
                    ++pos_;
                while (pos_ < text_.size() && std::isdigit(static_cast<unsigned char>(text_[pos_])))
                    ++pos_;
            }

            return std::stod(text_.substr(start, pos_ - start));
        }

        bool matchLiteral(const char *literal)
        {
            size_t length = std::char_traits<char>::length(literal);
            if (text_.compare(pos_, length, literal) == 0)
            {
                pos_ += length;
                return true;
            }
            return false;
        }

        void skipWhitespace()
        {
            while (pos_ < text_.size() && std::isspace(static_cast<unsigned char>(text_[pos_])))
                ++pos_;
        }

        bool peek(char ch) const
        {
            return pos_ < text_.size() && text_[pos_] == ch;
        }

        void expect(char ch)
        {
            skipWhitespace();
            if (pos_ >= text_.size() || text_[pos_] != ch)
                throw std::runtime_error("Unexpected JSON syntax.");
            ++pos_;
        }

        std::string text_;
        size_t pos_;
    };

    const JsonValue &requireField(const JsonObject &object, const std::string &name)
    {
        auto it = object.find(name);
        if (it == object.end())
            throw std::runtime_error("Missing required JSON field: " + name);
        return it->second;
    }

    glm::dvec3 parseVec3(const JsonValue &value)
    {
        const JsonArray &array = value.asArray();
        if (array.size() != 3)
            throw std::runtime_error("Control point must contain exactly 3 coordinates.");
        return glm::dvec3(array[0].asNumber(), array[1].asNumber(), array[2].asNumber());
    }

    std::vector<double> parseNumberArray(const JsonValue &value)
    {
        const JsonArray &array = value.asArray();
        std::vector<double> result;
        result.reserve(array.size());
        for (const JsonValue &entry : array)
            result.push_back(entry.asNumber());
        return result;
    }

    std::vector<glm::dvec3> parseControlPoints(const JsonValue &value)
    {
        const JsonArray &array = value.asArray();
        std::vector<glm::dvec3> result;
        result.reserve(array.size());
        for (const JsonValue &entry : array)
            result.push_back(parseVec3(entry));
        return result;
    }

    tinynurbs::RationalCurve3d parseCurve(const JsonObject &object)
    {
        tinynurbs::RationalCurve3d curve;
        curve.degree = static_cast<unsigned int>(requireField(object, "degree").asNumber());
        curve.knots = parseNumberArray(requireField(object, "knots"));
        curve.control_points = parseControlPoints(requireField(object, "control_points"));
        curve.weights = parseNumberArray(requireField(object, "weights"));

        if (curve.control_points.size() != curve.weights.size())
            throw std::runtime_error("Curve control point count must match weights count.");
        if (curve.knots.size() != curve.control_points.size() + curve.degree + 1)
            throw std::runtime_error("Curve knot vector size is inconsistent with degree/control points.");

        return curve;
    }
}

TestCaseLibrary loadTestCaseLibrary(const std::filesystem::path &jsonPath)
{
    std::ifstream input(jsonPath);
    if (!input)
        throw std::runtime_error("Failed to open test case file: " + jsonPath.string());

    std::ostringstream buffer;
    buffer << input.rdbuf();

    JsonParser parser(buffer.str());
    JsonValue root = parser.parse();
    const JsonObject &rootObject = root.asObject();

    TestCaseLibrary library;
    library.defaultCase = requireField(rootObject, "default_case").asString();

    const JsonArray &cases = requireField(rootObject, "cases").asArray();
    library.cases.reserve(cases.size());

    for (const JsonValue &caseValue : cases)
    {
        const JsonObject &caseObject = caseValue.asObject();
        TestCaseData data;
        data.name = requireField(caseObject, "name").asString();

        auto notesIt = caseObject.find("notes");
        if (notesIt != caseObject.end())
            data.notes = notesIt->second.asString();

        const JsonArray &curveArray = requireField(caseObject, "curves").asArray();
        data.curves.reserve(curveArray.size());
        for (const JsonValue &curveValue : curveArray)
            data.curves.push_back(parseCurve(curveValue.asObject()));

        if (data.curves.empty())
            throw std::runtime_error("Each test case must contain at least one curve.");

        library.cases.push_back(std::move(data));
    }

    return library;
}

const TestCaseData &findTestCase(const TestCaseLibrary &library, const std::string &caseName)
{
    for (const TestCaseData &data : library.cases)
    {
        if (data.name == caseName)
            return data;
    }
    throw std::runtime_error("Test case not found: " + caseName);
}
