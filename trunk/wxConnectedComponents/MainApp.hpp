/*!
 * \file
 * \author Chris de Pujo
 */

#ifndef __MAINAPP__
#define __MAINAPP__

#include <wx/app.h>

#include "MainFrame.hpp"
#include "MainPanel.hpp"

class MainApp : public wxApp
{
public:

    virtual bool OnInit();

    void OnIdle(wxIdleEvent & event);

    MainFrame * frame;
    MainPanel * panel;

    DECLARE_EVENT_TABLE()
};

#endif // __MAINAPP__
