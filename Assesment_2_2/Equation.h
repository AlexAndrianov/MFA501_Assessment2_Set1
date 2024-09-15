#pragma once
#include <functional>
#include <memory>
#include <unordered_map>
#include <string>
#include <stdexcept>
#include <list>
#include <stack>
#include <deque>
#include <optional>
#include <math.h>

namespace math {

enum class TokenType : int {
    sin = 0,
    cos,
    tan,
    cat,
    asin,
    acos,
    atan,
    actan,
    exp,
    ln,
    log,
    ext,
    sqrt,
    plus,
    minus,
    multipl,
    devision,
    left_bracket,
    right_bracket,
    value,
    var,
    module,
    all,
    any,
    bracket_gr,
    module_gr,
    sin_gr,
    cos_gr,
    tan_gr,
    cat_gr,
    asin_gr,
    acos_gr,
    atan_gr,
    actan_gr,
    exp_gr,
    ln_gr,
    log_gr,
    ext_gr,
    sqrt_gr,
    plus_gr,
    minus_gr,
    multipl_gr,
    devision_gr,
    single_minus_gr,
    nothing,
};

const static std::unordered_map<std::string, TokenType> g_TokenLiterals = {
    {"sin", TokenType::sin },
    {"cos", TokenType::cos },
    {"tan", TokenType::tan },
    {"cat", TokenType::cat },
    {"asin", TokenType::asin },
    {"acos", TokenType::acos },
    {"atan", TokenType::atan },
    {"actan", TokenType::actan },
    {"exp", TokenType::exp },
    {"ln", TokenType::ln },
    {"log", TokenType::log },
    {"^", TokenType::ext },
    {"sqrt", TokenType::sqrt },
    {"+", TokenType::plus },
    {"-", TokenType::minus },
    {"*", TokenType::multipl },
    {"/", TokenType::devision },
    {"(", TokenType::left_bracket },
    {")", TokenType::right_bracket },
    {"x", TokenType::var },
    {"|", TokenType::module }
};

const static std::unordered_map<TokenType, std::string> g_LiteralTokens = [](){
    std::unordered_map<TokenType, std::string> res;

    for(const auto &tokenToLiteral : g_TokenLiterals)
        res[tokenToLiteral.second] = tokenToLiteral.first;

    return res;
}();

class CalculationContext {
public:
    CalculationContext(double parameterV)
        : _parameterV(parameterV) {}

    virtual ~CalculationContext() = default;

    double _parameterV = 0;
};

class Operator {
public:
    virtual ~Operator() = default;

    virtual double produce(const CalculationContext &context) const = 0;
    virtual std::string toString() const = 0;
};

class ConstantOperator : public Operator{
public:

    ConstantOperator(double v) : _v(v) {}

    virtual double produce(const CalculationContext &) const
    {
        return _v;
    }

    virtual std::string toString() const
    {
        return std::to_string(_v);
    }

    double _v = 0;
};

class VariableOperator : public Operator{
public:

    virtual double produce(const CalculationContext &context) const
    {
        return context._parameterV;
    }

    virtual std::string toString() const
    {
        return "x";
    }
};

class BinaryOperator : public Operator{
public:

    BinaryOperator(TokenType t,
                   std::shared_ptr<Operator> left,
                   std::shared_ptr<Operator> right,
                   std::function<double(double lV, double rV)> action) :
        _t(t), _left(left), _right(right), _action(action)
    { }

    virtual double produce(const CalculationContext &context) const
    {
        const auto lVl = _left->produce(context);
        const auto rVl = _right->produce(context);

        return _action(lVl, rVl);
    }

    virtual std::string toString() const
    {
        if(_t == TokenType::log)
            return g_LiteralTokens.at(_t) + _left->toString() + _right->toString();

        return _left->toString() + g_LiteralTokens.at(_t) + _right->toString();
    }

    TokenType _t;
    std::shared_ptr<Operator> _left;
    std::shared_ptr<Operator> _right;
    std::function<double(double lV, double rV)> _action;
};

class UnaryOperator : public Operator {
public:

    UnaryOperator(TokenType t,
                  std::shared_ptr<Operator> sub_group,
                  std::function<double(double v)> action)
        : _t(t), _sub_group(sub_group), _action(action)
    {}

    virtual double produce(const CalculationContext &context) const
    {
        return _action(_sub_group->produce(context));
    }

    virtual std::string toString() const
    {
        switch(_t)
        {
        case TokenType::bracket_gr:
            return "(" + _sub_group->toString() + ")";
        case TokenType::module_gr:
            return "|" + _sub_group->toString() + "|";
        default:
            return g_LiteralTokens.at(_t) + _sub_group->toString();
        }
    }

    TokenType _t;
    std::shared_ptr<Operator> _sub_group;
    std::function<double(double v)> _action;
};

struct Token {
    Token(TokenType type) : _type(type) {}
    virtual ~Token() = default;

    TokenType _type;
};

struct TokenValue : public Token {
    TokenValue(TokenType type, double v) : Token(type), _v(v) {}

    double _v = 0.;
};

struct TokenGroup : public Token {
    TokenGroup(TokenType type) : Token(type) {}

    std::list<std::shared_ptr<Token>> _group;
};

class Equation {
    std::list<std::shared_ptr<Token>> tokenize(const std::string &script) const
    {
        std::list<std::shared_ptr<Token>> res;

        bool isValue = false;
        std::string buffer;

        for(auto ch : script)
        {
            if(std::isdigit(ch) || ch == '.' || ch == ',')
            {
                if(!buffer.empty() && !isValue)
                {
                    throw std::runtime_error("Undefined symbol: " + buffer);
                }

                isValue = true;
            }
            else if(isValue)
            {
                res.emplace_back(std::make_shared<TokenValue>(TokenType::value, std::stod(buffer)));
                isValue = false;
                buffer.clear();
            }

            buffer.push_back(ch);

            auto it = g_TokenLiterals.find(buffer);
            if(it != g_TokenLiterals.end())
            {
                res.emplace_back(std::make_shared<Token>(it->second));
                buffer.clear();
            }
        }

        if(isValue)
        {
            res.emplace_back(std::make_shared<TokenValue>(TokenType::value, std::stod(buffer)));
        }
        else if(!buffer.empty())
            throw std::runtime_error("Undefined symbol: " + buffer);

        return res;
    }

    std::list<std::shared_ptr<Token>> grammatics_rules_apply(std::list<std::shared_ptr<Token>> tokens)
    {
        struct Grammatic {
            TokenType _grType;
            std::deque<TokenType> _tokenOrdering;
        };

        // order is valuable
        static const std::list<Grammatic> s_Grammatics = {
            { TokenType::bracket_gr, {TokenType::left_bracket, TokenType::all, TokenType::right_bracket} },
            { TokenType::module_gr, {TokenType::module, TokenType::all, TokenType::module} },
            { TokenType::asin_gr, {TokenType::asin, TokenType::bracket_gr}},
            { TokenType::acos_gr, {TokenType::acos, TokenType::bracket_gr}},
            { TokenType::atan_gr, {TokenType::atan, TokenType::bracket_gr}},
            { TokenType::sin_gr, {TokenType::sin, TokenType::bracket_gr}},
            { TokenType::cos_gr, {TokenType::cos, TokenType::bracket_gr}},
            { TokenType::tan_gr, {TokenType::tan, TokenType::bracket_gr}},
            { TokenType::cat_gr, {TokenType::cat, TokenType::bracket_gr}},
            { TokenType::exp_gr, {TokenType::exp, TokenType::bracket_gr}},
            { TokenType::ln_gr, {TokenType::ln, TokenType::bracket_gr}},
            { TokenType::sqrt_gr, {TokenType::sqrt, TokenType::bracket_gr}},
            { TokenType::single_minus_gr, {TokenType::nothing, TokenType::minus, TokenType::any}},
            { TokenType::single_minus_gr, {TokenType::plus, TokenType::minus, TokenType::any}},
            { TokenType::single_minus_gr, {TokenType::multipl, TokenType::minus, TokenType::any}},
            { TokenType::single_minus_gr, {TokenType::minus, TokenType::minus, TokenType::any}},
            { TokenType::single_minus_gr, {TokenType::devision, TokenType::minus, TokenType::any}},
            { TokenType::single_minus_gr, {TokenType::ext, TokenType::minus, TokenType::any}},
            { TokenType::single_minus_gr, {TokenType::log, TokenType::minus, TokenType::any}},
            { TokenType::ext_gr, {TokenType::any, TokenType::ext, TokenType::any}},
            { TokenType::log_gr, {TokenType::log, TokenType::any, TokenType::bracket_gr}},
            { TokenType::multipl_gr, {TokenType::any, TokenType::multipl, TokenType::any}},
            { TokenType::devision_gr, {TokenType::any, TokenType::devision, TokenType::any}},
            { TokenType::plus_gr, {TokenType::any, TokenType::plus, TokenType::any}},
            { TokenType::minus_gr, {TokenType::any, TokenType::minus, TokenType::any}}
        };

        auto ruleProcessor = [this, &tokens](const Grammatic &rule){
            std::stack<TokenType> expectedTokens;
            TokenType expectedToken;

            decltype(tokens.end()) firstSubToken;
            decltype(tokens.end()) endSubToken;

            auto resetExpectation = [&](){
                expectedTokens = std::stack<TokenType>(rule._tokenOrdering);
                expectedToken = expectedTokens.top();
                expectedTokens.pop();

                firstSubToken = tokens.end();
                endSubToken = tokens.end();
            };

            resetExpectation();

            bool allMode = false;
            auto sub_bracket_index = 0;
            TokenType bracket_type;

            for(auto rit = tokens.rbegin(); rit != tokens.rend(); rit++)
            {
                const auto currentNodeType = (*rit)->_type;

                if(expectedToken != TokenType::any && currentNodeType != expectedToken)
                {
                    if(!allMode && rule._tokenOrdering.size() != expectedTokens.size() + 1)
                        resetExpectation();

                    if(allMode && currentNodeType == bracket_type)
                        sub_bracket_index++;

                    continue;
                }

                if(expectedTokens.empty())
                {
                    if(allMode && sub_bracket_index)
                    {
                        sub_bracket_index--;
                        continue;
                    }

                    firstSubToken = rit.base();
                    --firstSubToken;

                    if(rule._grType == TokenType::single_minus_gr)
                        firstSubToken++;

                    break;
                }

                if(rule._tokenOrdering.size() == expectedTokens.size() + 1)
                {
                    endSubToken = rit.base();
                }

                auto previousToken = expectedToken;
                expectedToken = expectedTokens.top();
                expectedTokens.pop();

                if(expectedToken == TokenType::all)
                {
                    bracket_type = previousToken;

                    expectedToken = expectedTokens.top();
                    expectedTokens.pop();

                    allMode = true;
                    sub_bracket_index = 0;
                }
            }

            if(expectedToken == TokenType::nothing)
                firstSubToken = tokens.begin();

            if(firstSubToken == tokens.end())
                return false;

            auto group = std::make_shared<TokenGroup>(rule._grType);

            if(allMode)
            {
                auto firstFactSubToken = firstSubToken;
                auto endFactSubToken = endSubToken;

                group->_group = grammatics_rules_apply(std::list<std::shared_ptr<Token>>{++firstFactSubToken, --endFactSubToken});
            }
            else
            {
                group->_group = std::list<std::shared_ptr<Token>>{firstSubToken, endSubToken};
            }

            tokens.erase(firstSubToken, endSubToken);
            tokens.insert(endSubToken, group);
            return true;
        };

        for(const auto &rule : s_Grammatics)
        {
            while(ruleProcessor(rule));
        }

        return tokens;
    }

    std::shared_ptr<Operator> tokenToOperator(std::shared_ptr<Token> token)
    {
        if(!token)
            throw std::runtime_error("Parser error.");

        if(token->_type == TokenType::var)
            return std::make_shared<VariableOperator>();

        if(auto tokenV = std::dynamic_pointer_cast<TokenValue>(token))
            return std::make_shared<ConstantOperator>(tokenV->_v);

        if(auto tokenG = std::dynamic_pointer_cast<TokenGroup>(token))
        {
            if(tokenG->_group.empty())
                return nullptr;

            switch (tokenG->_type) {
            case TokenType::bracket_gr:
                return std::make_shared<UnaryOperator>(tokenG->_type, tokenToOperator(*tokenG->_group.begin()),
                                                       [](double v){ return v; });
            case TokenType::module_gr:
                return std::make_shared<UnaryOperator>(tokenG->_type, tokenToOperator(*tokenG->_group.begin()),
                                                       [](double v){ return v < 0 ? -v : v; });
            case TokenType::asin_gr:
                return std::make_shared<UnaryOperator>((*tokenG->_group.begin())->_type, tokenToOperator(*tokenG->_group.rbegin()),
                                                       [](double v){ return asin(v); });
            case TokenType::acos_gr:
                return std::make_shared<UnaryOperator>((*tokenG->_group.begin())->_type, tokenToOperator(*tokenG->_group.rbegin()),
                                                       [](double v){ return acos(v); });
            case TokenType::atan_gr:
                return std::make_shared<UnaryOperator>((*tokenG->_group.begin())->_type, tokenToOperator(*tokenG->_group.rbegin()),
                                                       [](double v){ return atan(v); });
            case TokenType::sin_gr:
                return std::make_shared<UnaryOperator>((*tokenG->_group.begin())->_type, tokenToOperator(*tokenG->_group.rbegin()),
                                                       [](double v){ return sin(v); });
            case TokenType::cos_gr:
                return std::make_shared<UnaryOperator>((*tokenG->_group.begin())->_type, tokenToOperator(*tokenG->_group.rbegin()),
                                                       [](double v){ return cos(v); });
            case TokenType::tan_gr:
                return std::make_shared<UnaryOperator>((*tokenG->_group.begin())->_type, tokenToOperator(*tokenG->_group.rbegin()),
                                                       [](double v){ return tan(v); });
            case TokenType::cat_gr:
                return std::make_shared<UnaryOperator>((*tokenG->_group.begin())->_type, tokenToOperator(*tokenG->_group.rbegin()),
                                                       [](double v){ return 1.0/tan(v); });
            case TokenType::exp_gr:
                return std::make_shared<UnaryOperator>((*tokenG->_group.begin())->_type, tokenToOperator(*tokenG->_group.rbegin()),
                                                       [](double v){ return exp(v); });
            case TokenType::ln_gr:
                return std::make_shared<UnaryOperator>((*tokenG->_group.begin())->_type, tokenToOperator(*tokenG->_group.rbegin()),
                                                       [](double v){ return log(v); });
            case TokenType::sqrt_gr:
                return std::make_shared<UnaryOperator>((*tokenG->_group.begin())->_type, tokenToOperator(*tokenG->_group.rbegin()),
                                                       [](double v){ return sqrt(v); });
            case TokenType::single_minus_gr:
                return std::make_shared<UnaryOperator>(TokenType::minus, tokenToOperator(*tokenG->_group.rbegin()),
                                                       [](double v){ return -v; });
            case TokenType::ext_gr:
                return std::make_shared<BinaryOperator>(TokenType::ext, tokenToOperator(*tokenG->_group.begin()),
                                                        tokenToOperator(*tokenG->_group.rbegin()),
                                                        [](double lv, double rv){ return pow(lv, rv); });
            case TokenType::log_gr:
                return std::make_shared<BinaryOperator>(TokenType::log, tokenToOperator(*(++(tokenG->_group.begin()))),
                                                        tokenToOperator(*tokenG->_group.rbegin()),
                                                        [](double lv, double rv){ return log2(rv) / log2(lv); });
            case TokenType::multipl_gr:
                return std::make_shared<BinaryOperator>((*(++(tokenG->_group.begin())))->_type, tokenToOperator(*tokenG->_group.begin()),
                                                        tokenToOperator(*tokenG->_group.rbegin()),
                                                        [](double lv, double rv){ return lv * rv; });
            case TokenType::devision_gr:
                return std::make_shared<BinaryOperator>((*(++(tokenG->_group.begin())))->_type, tokenToOperator(*tokenG->_group.begin()),
                                                        tokenToOperator(*tokenG->_group.rbegin()),
                                                        [](double lv, double rv){ return lv/rv; });
            case TokenType::plus_gr:
                return std::make_shared<BinaryOperator>((*(++(tokenG->_group.begin())))->_type, tokenToOperator(*tokenG->_group.begin()),
                                                        tokenToOperator(*tokenG->_group.rbegin()),
                                                        [](double lv, double rv){ return lv + rv; });
            case TokenType::minus_gr:
                return std::make_shared<BinaryOperator>((*(++(tokenG->_group.begin())))->_type, tokenToOperator(*tokenG->_group.begin()),
                                                        tokenToOperator(*tokenG->_group.rbegin()),
                                                        [](double lv, double rv){ return lv - rv; });
            default:
                throw std::runtime_error("Parser error.");
            }
        }

        throw std::runtime_error("Undefined token: " + std::to_string(static_cast<int>(token->_type)));
    }

public:

    // -logx(77)*sin(x^(ln(x+atan((sqrt(x)+-44.4645651))/x*exp(2*x))))
    void parse(const std::string &script)
    {
        auto tokens = tokenize(script);
        tokens = grammatics_rules_apply(tokens);

        if(tokens.size() != 1)
            throw std::runtime_error("Mistaken equation: " + script);

        _sintaxis_tree_root = tokenToOperator(*tokens.begin());
    }

    double calculateIntegral(double from, double to, double tol = 1e-6)
    {
        const auto startDevision = 10;
        const auto multiplicator = 10;
        auto iterNumber = startDevision;

        std::optional<double> previousResult = std::nullopt;
        double res = 0;

        while(true)
        {
            const float dx = fabs(to - from) / iterNumber;

            auto x = from;
            auto y1 = _sintaxis_tree_root->produce({x});
            res = 0;

            for(auto i = 0; i < iterNumber; i++)
            {
                x += dx;
                auto y2 = _sintaxis_tree_root->produce({x});

                const auto y = (y1 + y2)/2.;
                y1 = y2;
                res += y*dx;
            }

            if(previousResult.has_value() && fabs(*previousResult - res) < tol)
                break;

            iterNumber *= multiplicator;
            previousResult = res;
        }

        return res;
    }

    std::shared_ptr<Operator> _sintaxis_tree_root;
};

}
