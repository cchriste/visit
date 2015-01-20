/*    
 * Copyright (C) 1999-2011 University of Connecticut Health Center    
 *    
 * Licensed under the MIT License (the "License").    
 * You may not use this file except in compliance with the License.    
 * You may obtain a copy of the License at:    
 *    
 *  http://www.opensource.org/licenses/mit-license.php    
 */    
    
#include "DivideByZeroException.h"    
    
using namespace VCell;

DivideByZeroException::DivideByZeroException(std::string msg) : ExpressionException("DivideByZeroException", msg)    
{    
}    
    
DivideByZeroException::~DivideByZeroException(void)    
{    
}    
