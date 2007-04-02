/*****************************************************************************
*
* Copyright (c) 2000 - 2006, The Regents of the University of California
* Produced at the Lawrence Livermore National Laboratory
* All rights reserved.
*
* This file is part of VisIt. For details, see http://www.llnl.gov/visit/. The
* full copyright notice is contained in the file COPYRIGHT located at the root
* of the VisIt distribution or at http://www.llnl.gov/visit/copyright.html.
*
* Redistribution  and  use  in  source  and  binary  forms,  with  or  without
* modification, are permitted provided that the following conditions are met:
*
*  - Redistributions of  source code must  retain the above  copyright notice,
*    this list of conditions and the disclaimer below.
*  - Redistributions in binary form must reproduce the above copyright notice,
*    this  list of  conditions  and  the  disclaimer (as noted below)  in  the
*    documentation and/or materials provided with the distribution.
*  - Neither the name of the UC/LLNL nor  the names of its contributors may be
*    used to  endorse or  promote products derived from  this software without
*    specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
* ARE  DISCLAIMED.  IN  NO  EVENT  SHALL  THE  REGENTS  OF  THE  UNIVERSITY OF
* CALIFORNIA, THE U.S.  DEPARTMENT  OF  ENERGY OR CONTRIBUTORS BE  LIABLE  FOR
* ANY  DIRECT,  INDIRECT,  INCIDENTAL,  SPECIAL,  EXEMPLARY,  OR CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT  LIMITED TO, PROCUREMENT OF  SUBSTITUTE GOODS OR
* SERVICES; LOSS OF  USE, DATA, OR PROFITS; OR  BUSINESS INTERRUPTION) HOWEVER
* CAUSED  AND  ON  ANY  THEORY  OF  LIABILITY,  WHETHER  IN  CONTRACT,  STRICT
* LIABILITY, OR TORT  (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY  WAY
* OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
* DAMAGE.
*
*****************************************************************************/

#include <VisItControlInterface_V1.h>
#include <VisItDataInterface_V1.h>
#include <SimHelperFunc.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* set to maximum length of command name string  */
#define MAX_CMD_STR_LEN 64

/* current number of active channels to the custom
   UI update information in the meta data. */
int VisIt_CurrentIndex = 0;
const int true = 1;
const int false = 0;

/*****************************************************************************
** Function: int VisItFindCMD ( char *name )
**
** Purpose:
**   Searches for a match between name and a name in the metat data custom
**   command array. If it does not find it findCMD returns -1.
**   If it does find a match findCMD return the index for the command
**
** Arguments:
**   mdd       : meta data structure holding the commands arrays
**   name      : name of the ui component to search for.
**   customCMD : is this command on the custom Ui or the main Sim window
**
** Programmer: Shelly Prevost
** Creation:   Tue Dec 6 2005
**
** Modifications:
**
** *****************************************************************************/

int VisItFindCMD (VisIt_SimulationMetaData mdd, char *name, int customCMD)
{
    int found = -1;
    int i;

    if ( customCMD )
    {
         for ( i =0; i < mdd.numCustomCommands; i++)
         {
            if ( !strcmp (mdd.customCommands[i].name, name))
            {
                found = i;
            }
         }
    }
    else
    {
         for ( i =0; i < mdd.numGenericCommands; i++)
         {
            if ( !strcmp (mdd.genericCommands[i].name, name))
            {

                found = i;
            }
         }
    }

    return found;
}

/*****************************************************************************
** Function: void VisItUpdateMainSimWinGUI ( VisIt_SimulationMetaData mdd, char *name,
**                                           char *data, int enable )
**
** Purpose:
**   Searches for a match between name and a name in the metat data GENERIC
**   command array. If it does not find no data is set.
**   If it does find a match the text field is updated
**
** Arguments:
**   mdd  : meta data structure holding the commands arrays
**   name : name of the ui component to search for.
**
** Programmer: Shelly Prevost
** Creation:   Wed. Aug 30 2006
**
** Modifications:
**
** **************************************************************************** */

void VisItUpdateMainSimWinGUI ( VisIt_SimulationMetaData mdd, char *name, char *data, int enable )
{
    int index;
    index = VisItFindCMD (mdd, name,false);

    if ( (index >= 0) && ( index < mdd.numGenericCommands))
    {
      if ( enable != -1 )
         mdd.genericCommands[index].enabled = enable;
      if ( data )
      {
         strncpy(mdd.genericCommands[index].text, data,  MAX_CMD_STR_LEN-1);
      }
    }
}

/*****************************************************************************
** Function: void VisItCreateCMD ( VisIt_SimulationMetaData mdd, char *name )
**
** Purpose:
**   Creates a entry in the meta data to allow updating the ui
**   component that is named name. The entry will not be created
**   if it exceeds the max number of commands or if it is already
**   in the array.
**
** Arguments:
**   mdd  : meta data structure holding the commands arrays
**   name : name of the ui component to search for.
**
** Programmer: Shelly Prevost
** Creation:   Tue Dec 6 2005
**
** Modifications:
**
** *****************************************************************************/
void VisItCreateCMD ( VisIt_SimulationMetaData mdd, char *name )
{
  int index;

  if ( VisIt_CurrentIndex < mdd.numCustomCommands)
  {
      index = VisItFindCMD (mdd, name,true );
      if ( index < 0 )
          strncpy(mdd.customCommands[VisIt_CurrentIndex].name, name,MAX_CMD_STR_LEN-1);
      VisIt_CurrentIndex++;
  }

}

/*****************************************************************************
** Function: void VisItInitCMD ( VisItSimulationControlCommand *cmd )
**
** Purpose:
**   This function initialize the meta data structure so that
**   it is in a know state.
**
** Arguments:
**   cmd  : the structure in need of initialization
**
** Programmer: Shelly Prevost
** Creation:   Tue Dec 6 2005
**
** Modifications:
**
** *****************************************************************************/
void VisItInitCMD ( VisIt_SimulationControlCommand *cmd )
{
    cmd->name = malloc(MAX_CMD_STR_LEN);
    strncpy ( cmd->name, "VISIT_NONE",MAX_CMD_STR_LEN);
    cmd->argType = VISIT_CMDARG_NONE;
    cmd->enabled = 1;
    cmd->signal= malloc(MAX_CMD_STR_LEN);
    strncpy ( cmd->name, "VISIT_NONE",MAX_CMD_STR_LEN);
    cmd->text= malloc(MAX_CMD_STR_LEN);
    strncpy ( cmd->text, "VISIT_NONE",MAX_CMD_STR_LEN);
    cmd->className = malloc(MAX_CMD_STR_LEN);
    cmd->className[0] = '\0';
    cmd->value = malloc(MAX_CMD_STR_LEN);
    strncpy ( cmd->value, "10",MAX_CMD_STR_LEN);
    cmd->uiType = malloc(MAX_CMD_STR_LEN);
    cmd->uiType[0] = '\0';
    cmd->isOn = 1;
    cmd->parent = malloc(MAX_CMD_STR_LEN);
    cmd->parent[0] = '\0';

}

/*****************************************************************************
** Function: void VisItSetCMDEnable (VisIt_SimulationMetaData mdd, char *name, int enabledCMD )
**
** Purpose:
**   this is a helper function to look up the ui component by
**   by name and then set the enable variable.
**
** Arguments:
**   mdd  : meta data structure holding the commands arrays
**   name   : name of the ui component to search for.
**   enable : 0 for false and 1 for true.
**
** Programmer: Shelly Prevost
** Creation:   Tue Dec 6 2005
**
** Modifications:
**
** **************************************************************************** */
void VisItSetCMDEnable (VisIt_SimulationMetaData mdd, char *name, int enabledCMD )
{
    int index;

    index = VisItFindCMD (mdd, name,true);
    if ( (index >= 0) && ( index < mdd.numCustomCommands))
    {
      mdd.customCommands[index].enabled = enabledCMD;
    }
}

/*****************************************************************************
** Function: void VisItSetCMDIsOn (VisIt_SimulationMetaData mdd, char *name, int isOn)
**
** Purpose:
**   this is a helper function to look up the ui component by
**   by name and then set the isOn variable.
**
** Arguments:
**   mdd  : meta data structure holding the commands arrays
**   name   : name of the ui component to search for.
**   enable : 0 for false and 1 for true.
**
** Programmer: Shelly Prevost
** Creation:   Tue Dec 6 2005
**
** Modifications:
**
** *****************************************************************************/
void VisItSetCMDIsOn (VisIt_SimulationMetaData mdd, char *name, int isOn )
{
    int index;

    index = VisItFindCMD ( mdd, name,true);
    if ( (index >= 0) && ( index < mdd.numCustomCommands))
    {
        mdd.customCommands[index].isOn = isOn;
    }
}

/*****************************************************************************
** Function: void VisItSetCMDValue (VisIt_SimulationMetaData mdd, char *name, int value )
**
** Purpose:
**   Searches for a match between name and a name in the metat data custom
**   command array. If it does not find a match it sets the meta data for
**   that entry to value.
**
** Arguments:
**   mdd  : meta data structure holding the commands arrays
**   name   : name of the ui component to search for.
**   value  : the value to set the found component too.
**
** Programmer: Shelly Prevost
** Creation:   Tue Dec 6 2005
**
** Modifications:
**
** *****************************************************************************/
void VisItSetCMDValue (VisIt_SimulationMetaData mdd, char *name, int value )
{
    int index;
    char strValue[MAX_CMD_STR_LEN];

    sprintf (  strValue, "%5d", value);
    index = VisItFindCMD ( mdd, name,true);
    if ( (index >= 0) && ( index < mdd.numCustomCommands))
    {
      strncpy(mdd.customCommands[index].value, strValue,  MAX_CMD_STR_LEN-1);
    }
}

/*****************************************************************************
** Function: void VisItSetCMDText (VisIt_SimulationMetaData mdd, char *name, char *text )
**
** Purpose:
**   Searches for a match between name and a name in the metat data custom
**   command array. If it does not find a match it sets the text value for
**   that ui entry too text.
**
** Arguments:
**   mdd  : meta data structure holding the commands arrays
**   name : name of the ui component to search for.
**   text : the text to change it too.
**
** Programmer: Shelly Prevost
** Creation:   Tue Dec 6 2005
**
** Modifications:
**
** *****************************************************************************/
void VisItSetCMDText (VisIt_SimulationMetaData mdd, char *name, char *text )
{
    int index;

    index = VisItFindCMD (mdd, name,true);
    if ( (index >= 0) && ( index < mdd.numCustomCommands))
    {
      strncpy(mdd.customCommands[index].text, text,  MAX_CMD_STR_LEN-1);
    }
}

/*****************************************************************************
** Function: void VisItInitGenericCMD (VisIt_SimulationMetaData &mdd, int index,char *name, char
**                                     *text,char *value, int enable )
**
** Purpose:
**   Searches for a match between name and a name in the metadata GENERIC
**   command array. If it does not find a match it does not set any data
**
** Arguments:
**   mdd  : meta data structure holding the commands arrays
**   name : index is the location in the generic command array to add this data too.
**   text : the text to change it too.
**   value: the value of the UI
**   enable: enable/disable the ui component
**
** Programmer: Shelly Prevost
** Creation:   Wed. Aug. 30 2006
**
** Modifications:
**
** *****************************************************************************/
void VisItInitGenericCMD (VisIt_SimulationMetaData mdd, int index,char *name, char *text,char *value, int enable )
{
    if ( (index >= 0) && ( index < mdd.numGenericCommands))
    {
      strncpy(mdd.genericCommands[index].name, name,  MAX_CMD_STR_LEN-1);
      strncpy(mdd.genericCommands[index].text, text,  MAX_CMD_STR_LEN-1);
      strncpy(mdd.genericCommands[index].value, value,  MAX_CMD_STR_LEN-1);
      mdd.genericCommands[index].argType = VISIT_CMDARG_NONE;
      mdd.genericCommands[index].enabled = enable;
    }
}

/*****************************************************************************
** Function: vvoid VisItLabelGenericButton (VisIt_SimulationMetaData mdd, int button,
**                                           char *text,int enable )
**
** Purpose:
**   Sets the <button> label to <text> and enable or disable the button according to
**   the value of <enable>
**
** Arguments:
**   mdd: meta data structure
**   button: the number of the button label you want to change. It must be less than 5.
**   text : the text to change it too.
**   enable: set the ui to be enabled
**
** Programmer: Shelly Prevost
** Creation:   Wed. Aug. 30 2006
**
** Modifications:
**
** *****************************************************************************/
void VisItLabelGenericButton (VisIt_SimulationMetaData *mdd, int button, char *text,int enable )
{
    if ( (button >= 0) && ( button < mdd->numGenericCommands))
    {
      strncpy(mdd->genericCommands[button].name, text,  MAX_CMD_STR_LEN-1);
      strncpy(mdd->genericCommands[button].text, text,  MAX_CMD_STR_LEN-1);
      if ( enable != -1 ) mdd->genericCommands[button].enabled = enable;
    }
}


/*****************************************************************************
** Function: VisItsetCMD (VisIt_SimulationMetaData mdd, VisItSimulationControlCommand cmd  )
**
** Purpose:
**   Searches for a match between cmd.name and a name in the metat data custom
**   command array. If it does find a match it sets all the values in that
**   entry to all the values in cmd
**
** Arguments:
**   mdd  : meta data structure holding the commands arrays
**   cmd : data structure to copy to the matching entry.
**
** Programmer: Shelly Prevost
** Creation:   Tue Dec 6 2005
**
** Modifications:
**
** *****************************************************************************/
void VisItSetCMD ( VisIt_SimulationMetaData mdd, VisIt_SimulationControlCommand cmd  )
{
    int index;

    index = VisItFindCMD ( mdd, cmd.name,true);
    if ( (index >= 0) && ( index < mdd.numCustomCommands))
    {
        mdd.customCommands[index].argType = cmd.argType;
        mdd.customCommands[index].enabled = cmd.enabled;
        mdd.customCommands[index].isOn = cmd.isOn;
        strncpy(mdd.customCommands[index].name,       cmd.name,     MAX_CMD_STR_LEN-1);
        strncpy(mdd.customCommands[index].signal,     cmd.signal,   MAX_CMD_STR_LEN-1);
        strncpy(mdd.customCommands[index].text,       cmd.text,     MAX_CMD_STR_LEN-1);
        strncpy(mdd.customCommands[index].className,  cmd.className,MAX_CMD_STR_LEN-1);
        strncpy(mdd.customCommands[index].value,      cmd.value,    MAX_CMD_STR_LEN-1);
        strncpy(mdd.customCommands[index].uiType,     cmd.uiType,   MAX_CMD_STR_LEN-1);
        strncpy(mdd.customCommands[index].parent,     cmd.parent,   MAX_CMD_STR_LEN-1);
    }
}



/*****************************************************************************
** Function: void VisItInitAllCMD(VisIt_SimulationMetaData *mdd, int MaxNumCustCMD  )
**
** Purpose:
**   This functions initializes the meta data so that all the command arrays are
**   initialized to some type of valid data. This prevent problems later on if
**   the meta is sent to VisIt before it has be set.
**
** Arguments:
**   mdd  : meta data structure holding the commands arrays
**   MaxNumCustCMD  : maximum number of command channels that will be open at
**   anyone time.
**
** Programmer: Shelly Prevost
** Creation:   Tue Dec 6 2005
**
** Modifications:
**
** **************************************************************************** */
void VisItInitAllCMD(VisIt_SimulationMetaData *mdd, int MaxNumCustCMD  )
{
   const int MaxNumGenericCommands = 20;
   int i;
   int index =0;
   /* These are the definitions for the           */
   /* generic commands found in the first window  */
   /* when VisIt open the first simulation window */
   mdd->numGenericCommands =MaxNumGenericCommands;
   mdd->genericCommands = malloc(sizeof(VisIt_SimulationControlCommand) *mdd->numGenericCommands);
   for ( i =0; i <mdd->numGenericCommands; i++) VisItInitCMD(&mdd->genericCommands[i]);

   VisItInitGenericCMD (*mdd, index++,"halt", "halt","NONE", 1 );
   VisItInitGenericCMD (*mdd, index++,"step", "step","NONE", 1 );
   VisItInitGenericCMD (*mdd, index++,"run", "run","NONE", 1 );
   VisItInitGenericCMD (*mdd, index++,"restart", "restart","NONE", 1 );
   VisItInitGenericCMD (*mdd, index++,"sync", "sync","NONE", 1 );
   VisItInitGenericCMD (*mdd, index++,"custom", "custom","NONE", 1 );
   VisItInitGenericCMD (*mdd, index++,"StartLabel", "Start","NONE", 1 );
   VisItInitGenericCMD (*mdd, index++,"StepLabel", "Step","NONE", 1 );
   VisItInitGenericCMD (*mdd, index++,"StopLabel", "Stop","NONE", 1 );
   VisItInitGenericCMD (*mdd, index++,"MessageViewerTextEdit", "-","NONE", 1 );
   VisItInitGenericCMD (*mdd, index++,"MessageViewerLabel", "Message Viewer","NONE", 1 );
   
   /* VisItInitGenericCMD (*mdd, index++,"StartLineEdit", "0.0000","10", 1 );
      VisItInitGenericCMD (*mdd, index++,"StepLineEdit", "1.0000","NONE", 1 );
      VisItInitGenericCMD (*mdd, index++,"StopLineEdit", "10.0000","NONE", 1 );  */
  
   mdd->numGenericCommands =index;
   mdd->numCustomCommands = MaxNumCustCMD;
   mdd->customCommands = malloc(sizeof(VisIt_SimulationControlCommand) *mdd->numCustomCommands);
 
   for ( i =0; i <mdd->numCustomCommands; i++) VisItInitCMD(&mdd->customCommands[i]);
}
