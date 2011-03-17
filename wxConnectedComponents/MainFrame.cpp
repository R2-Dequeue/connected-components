/*!
 * \file
 * \author Chris de Pujo
 */

#include "MainFrame.hpp"

#include <wx/wx.h>

BEGIN_EVENT_TABLE(MainFrame, wxFrame)
    EVT_MENU(MenuIDQuit, MainFrame::OnQuit)
    EVT_MENU(MenuIDAbout, MainFrame::OnAbout)
END_EVENT_TABLE()

MainFrame::MainFrame(const wxString & title,
                     const wxPoint & pos, // wxDefaultPosition
                     const wxSize & size,
                     long style)
    : wxFrame(NULL, wxID_ANY, title, pos, size, style)
{
    wxMenu * fileMenu = new wxMenu;

    fileMenu->Append(MenuIDAbout, _("&About"));
    fileMenu->AppendSeparator();
    fileMenu->Append(MenuIDQuit, _("E&xit"));

    wxMenuBar * menuBar = new wxMenuBar;
    menuBar->Append(fileMenu, _("&File"));

    SetMenuBar(menuBar);

    CreateStatusBar();
    SetStatusText(_("status text"));
}

void MainFrame::OnQuit(wxCommandEvent & WXUNUSED(event))
{
    Close(true);
}

void MainFrame::OnAbout(wxCommandEvent & WXUNUSED(event))
{
    wxMessageBox(_(""),
                 _("About wxConnectedComponents"),
                 wxOK | wxICON_INFORMATION,
                 this);
}
