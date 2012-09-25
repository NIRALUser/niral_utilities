/*******************************************************************
* FILE:     BatchControls.h
* PROJECT:  None
* AUTHORS:  None 
* DATE:     None 
* COMMENTS: None
*******************************************************************/
#include "BatchControls.h"

#include <qfiledialog.h>
#include <qmessagebox.h>
#include <qcombobox.h>
#include <qcheckbox.h>

BatchControls::BatchControls( QWidget* parent, const char* name, bool modal, WFlags fl)
 : Batch( parent, name, fl ) 
{
	m_rescaler = new ImageIntensityNormalizer();
}

BatchControls::~BatchControls()
{ 
}




void BatchControls::Rescale()
{
	double newmax;

	ImageIntensityNormalizer::ImagePointer m_image;
	ImageIntensityNormalizer::ImagePointer m_target_window;
	ImageIntensityNormalizer::ImagePointer m_source_window;
	ImageIntensityNormalizer::VectorList* m_vectorlist;
	
	ImageIntensityNormalizer::ImagePointer m_source;
	ImageIntensityNormalizer::ImagePointer m_sourceseg;
	ImageIntensityNormalizer::ImagePointer m_target;
	ImageIntensityNormalizer::ImagePointer m_targetseg;


	m_target = m_rescaler->ReadImage(gui_imgtarget->text());
	m_targetseg = m_rescaler->ReadImage(gui_segtarget->text());

	//Compute Mean values for target
	m_vectorlist = m_rescaler->ComputeMax(m_target,m_targetseg);

	//Compute new image;
	m_target_window = m_rescaler->TargetIntensityWindowing(m_target,*m_vectorlist,atof(g_sigma->text().latin1()),&newmax);

	//Save target 
	QString filename = gui_imgtarget->text().mid(gui_imgtarget->text().findRev("/")+1);
	QString dirname = g_outputdir->text() + filename.left(filename.findRev(".")) + g_suffix->text(); 
	m_rescaler->SaveImage(m_target_window,dirname.latin1());

	QListViewItem * item= gui_imgsource->firstChild();
    while( item ) 
	{
		std::cout << "Computing: " << item->text(0).latin1() << std::endl;
		m_source = m_rescaler->ReadImage(item->text(0).latin1());
		m_sourceseg = m_rescaler->ReadImage(item->text(1).latin1());

		//Compute Mean Values for source
		m_vectorlist = m_rescaler->ComputeMax(m_source,m_sourceseg);

		//Copmute new image for source
		m_source_window = m_rescaler->IntensityWindowing(m_source,*m_vectorlist,atof(g_sigma->text().latin1()),newmax);

		QString filename = item->text(0).mid(item->text(0).findRev("/")+1);
		QString dirname = g_outputdir->text() + filename.left(filename.findRev(".")) + g_suffix->text(); 
	
		if (g_adjust->isChecked())
		{
			m_image = m_rescaler->AdjustClasses(m_target_window,m_targetseg,m_source_window,m_sourceseg);	
			m_rescaler->SaveImage(m_image,dirname.latin1());
		}
		else
		{
			m_rescaler->SaveImage(m_source_window,dirname.latin1());
		}


		item = item->nextSibling();
	}

}

void BatchControls::SelectSource()
{
	QStringList filelist( QFileDialog::getOpenFileNames( "Images (*.gipl *.hdr *.v0 *.mha *.img *.raw *.gipl.gz *.hdr.gz *.v0.gz *.mha.gz *.img.gz *.raw.gz)",QString::null, this ) );
	
	filelist.sort();
	for ( QStringList::Iterator it = filelist.begin(); it != filelist.end(); ++it ) 
	{
		QString img = (*it).latin1();
		++it;
		if (it == filelist.end()) return;
		QString seg = (*it).latin1();
		QListViewItem * item = new QListViewItem(gui_imgsource, 0 );
		item->setText( 0, img);
	    item->setText( 1, seg);
		gui_imgsource->insertItem(item);
	}
}

void BatchControls::SelectTargetImage()
{
	QString m_open( QFileDialog::getOpenFileName( "*.gipl",QString::null, this ) );
	if (m_open.isEmpty())
		return;

	gui_imgtarget->setText(m_open);
}

void BatchControls::SelectTargetSegmentation()
{
	QString m_open( QFileDialog::getOpenFileName( QString::null,"*.gipl", this ) );
	if (m_open.isEmpty())
		return;

	gui_segtarget->setText(m_open);
}

void BatchControls::SelectOutputDir()
{
	QString m_outputdir( QFileDialog::getExistingDirectory( QString::null, this ) );
	if (m_outputdir.isEmpty())
		return;

	g_outputdir->setText(m_outputdir);
}


void BatchControls::AddLabel()
{
	if (!g_label->text().isEmpty())
	{
		m_rescaler->AddLabel(atoi(g_label->text().latin1()));
		g_labellist->insertItem(g_label->text());
	}
	else
	{
		m_rescaler->ListLabel(m_rescaler->ReadImage(gui_segtarget->text()));
	}
}


void BatchControls::RemoveLabel()
{
	
}

void BatchControls::Show()
{
	show();
}

