/*!
 * \file
 * \author Chris de Pujo
 */

#include "MainApp.hpp"

#include <wx/wx.h>

IMPLEMENT_APP(MainApp)

BEGIN_EVENT_TABLE(MainApp, wxApp)
    EVT_IDLE(MainApp::OnIdle)
END_EVENT_TABLE()

bool MainApp::OnInit()
{
    this->frame = new MainFrame(_("ConnectedComponents"),
                                wxPoint(50, 50),
                                wxSize(500, 500),
                                wxFULL_REPAINT_ON_RESIZE | wxDEFAULT_FRAME_STYLE);
    this->panel = new MainPanel(this->frame);
    wxTextCtrl * input1 = new wxTextCtrl(this->frame, wxID_ANY,
                                         _("x^2 + y^2 - 4"));
    wxTextCtrl * input2 = new wxTextCtrl(this->frame, wxID_ANY, _("x + y"));
    wxButton * button = new wxButton(this->frame, wxID_ANY, _("Graph"));

    wxBoxSizer * L1sizer = new wxBoxSizer(wxVERTICAL);
    wxBoxSizer * L2sizer = new wxBoxSizer(wxHORIZONTAL);
    wxBoxSizer * L3sizer = new wxBoxSizer(wxVERTICAL);

    L3sizer->Add(input1, 0, wxEXPAND | wxBOTTOM, 2);
    L3sizer->Add(input2, 0, wxEXPAND);

    L2sizer->Add(L3sizer, 1, wxRIGHT, 2);
    L2sizer->Add(button, 0, wxEXPAND | wxLEFT, 2);

    L1sizer->Add(this->panel, 1, wxEXPAND);
    L1sizer->Add(L2sizer, 0, wxEXPAND | wxALL, 2);

    frame->SetSizer(L1sizer);
    frame->SetAutoLayout(true);

    frame->Show(true);
    //SetTopWindow(frame); // Is this necessary?

    return true;
}

void MainApp::OnIdle(wxIdleEvent & event)
{
    if (panel != NULL)
        if (panel->moreLinesToDraw())
        {
            panel->paintNow();
            event.RequestMore();
        }
}
