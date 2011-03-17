/*!
 * \file
 * \author Chris de Pujo
 */

#ifndef __MAINFRAME__
#define __MAINFRAME__

#include <wx/wx.h>

class MainFrame : public wxFrame
{
public:

    MainFrame(const wxString & title,
              const wxPoint & pos,
              const wxSize & size,
              long style = wxDEFAULT_FRAME_STYLE);
    void OnQuit(wxCommandEvent & event);
    void OnAbout(wxCommandEvent & event);

    DECLARE_EVENT_TABLE()
};

enum
{
// idMenuAbout
// menuIDAbout
// menuIdAbout
// idAbout
// idAboutMenu
// MenuIDAbout
// MenuIdAbout
    MenuIDQuit = 1,
    MenuIDAbout
};

#endif // __MAINFRAME__
