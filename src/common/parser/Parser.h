#ifndef PARSER_H
#define PARSER_H
#include <parser_exports.h>

#include <vector>
#include <visitstream.h>

#include <Token.h>
#include <ParseTreeNode.h>
#include <Symbol.h>
#include <State.h>
#include <Grammar.h>
#include <TokenParseTreeNode.h>

// ****************************************************************************
//  Class:  ParseElem
//
//  Purpose:
//    An element of a parse stack.
//
//  Programmer:  Jeremy Meredith
//  Creation:    April  5, 2002
//
//  Modifications:
//    Jeremy Meredith, Wed Nov 24 11:59:11 PST 2004
//    Major refactoring.  Added storage of Token here because it no longer
//    derives directly from a parse tree node.
//
// ****************************************************************************
struct PARSER_API ParseElem
{
    const Symbol   *sym;
    ParseTreeNode  *node;
    Token          *token;

    ParseElem(Token * t)
    {
        sym = Symbol::Get(t->GetType());
        node = new TokenParseTreeNode(t);
        token = t;
    }
    ParseElem(const Symbol * s, ParseTreeNode * n)
    {
        sym = s;
        node = n;
        token = NULL;
    }
};

// ****************************************************************************
//  Class:  Parser
//
//  Purpose:
//    Abstract base for a parser.
//
//  Programmer:  Jeremy Meredith
//  Creation:    April  5, 2002
//
//  Modifications:
//    Jeremy Meredith, Wed Nov 24 11:59:55 PST 2004
//    Removed the ParserInterface base class (it wasn't doing much good).
//    Renamed this class Parser for consistency.  Added list of tokens
//    to the rule reduction method to make life easier for the implementor.
//
// ****************************************************************************
class PARSER_API Parser
{
public:
    Parser();
    virtual ~Parser() {}
    void    Init();
    void    ParseOneToken(Token *);
    bool    Accept() { return accept; }
    virtual ParseTreeNode *Parse(const std::string &) = 0;
    ParseTreeNode *GetParseTree() { return elems[0].node; }
    void    SetGrammar(Grammar * g) { G = g; }

protected:
    Grammar *G;
    std::vector<int> states;
    std::vector<ParseElem> elems;
    bool    accept;

protected:
    void    Shift(Token *, int);
    void    Reduce(int);
    virtual ParseTreeNode *ApplyRule(const Symbol &, const Rule *,
                                   std::vector<ParseTreeNode*>&,
                                   std::vector<Token*>&,
                                   Pos) = 0;
    void    PrintState(ostream &);
};

#endif
