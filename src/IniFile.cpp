#include "IniFile.h"
#include "Utilities.h"
#include "IO.h"

CIniFile::CIniFile() = default;

// A function to trim whitespace from both sides of a given string
void CIniFile::Trim(std::string& str, const std::string& CharsToTrim, int TrimDir)
{
    size_t startIndex = str.find_first_not_of(CharsToTrim);
    if (startIndex == std::string::npos){ str.erase(); return; }
    if (TrimDir < 2) { str = str.substr(startIndex, str.size()-startIndex); }
    if (TrimDir!=1) { str = str.substr(0, str.find_last_not_of(CharsToTrim) + 1); }
}

bool CIniFile::Load(const std::string& FileName, std::vector<Record>& content)
{
	std::string s;																// Holds the current line from the ini file
	std::string CurrentSection;													// Holds the current section name

	std::ifstream inFile (FileName.c_str());										// Create an input filestream
	if (!inFile.is_open()) { return false; }									// If the input file doesn't open, then return
	content.clear();														// Clear the content vector

	std::string comments;													// A string to store comments in

	while(!std::getline(inFile, s).eof())									// Read until the end of the file
	{
		Trim(s);															// Trim whitespace from the ends
		if(!s.empty())														// Make sure its not a blank line
		{
			Record r{};														// Define a new record

			if((s[0]=='#')||(s[0]==';'))									// Is this a commented line?
			{
				if ((s.find('[')==std::string::npos)&&							// If there is no [ or =
					(s.find('=')==std::string::npos))							// Then it's a comment
				{
					comments += s + '\n';									// Add the comment to the current comments string
				} else {
					r.Commented = s[0];										// Save the comment character
					s.erase(s.begin());										// Remove the comment for further processing
					Trim(s);
				}// Remove any more whitespace
			} else { r.Commented = ' '; }										// else mark it as not being a comment

			if(s.find('[')!=std::string::npos)									// Is this line a section?
			{		
				s.erase(s.begin());											// Erase the leading bracket
				s.erase(s.find(']'));										// Erase the trailing bracket
				r.Comments = comments;										// Add the comments string (if any)
				comments = "";												// Clear the comments for re-use
				r.Section = s;												// Set the Section value
				r.Key = "";													// Set the Key value
				r.Value = "";												// Set the Value value
				CurrentSection = s;
			}

			if(s.find('=')!=std::string::npos)									// Is this line a Key/Value?
			{
				r.Comments = comments;										// Add the comments string (if any)
				comments = "";												// Clear the comments for re-use
				r.Section = CurrentSection;									// Set the section to the current Section
				r.Key = s.substr(0,s.find('='));							// Set the Key value to everything before the = sign
				r.Value = s.substr(s.find('=')+1);							// Set the Value to everything after the = sign
			}
			if(comments.empty())	{											// Don't add a record yet if its a comment line
				content.push_back(r);										// Add the record to content
			}
		}
	}
	
	inFile.close();															// Close the file
	return true;
}

bool CIniFile::Save(const std::string& FileName, const std::vector<Record>& content)
{
	std::ofstream outFile (FileName.c_str());									// Create an output filestream
	if (!outFile.is_open()) { return false; }									// If the output file doesn't open, then return

	for (const auto& iContent : content)				// Loop through each vector
	{
		outFile << iContent.Comments;										// Write out the comments
		if(iContent.Key.empty())	{										// Is this a section?
			outFile << iContent.Commented << "[" 
			<< iContent.Section << "]" << std::endl;							// Then format the section
		}
		else {
			outFile << iContent.Commented << iContent.Key  
			<< "=" << iContent.Value << std::endl;								// Else format a key/value
		}
	}

	outFile.close();														// Close the file
	return true;
}

std::string CIniFile::Content(const std::string& FileName)
{
	std::string s;															// Hold our return string
	std::vector<Record> content;													// Holds the current record													// Holds the current record

	if (Load(FileName, content))											// Make sure the file loads
	{
		int i = 0;
		for (const auto& iContent : content)								// Loop through the content
		{
			if(!iContent.Comments.empty()) { s += iContent.Comments; }			// Add the comments
			if(iContent.Commented != ' ') { s += iContent.Commented; }		// If this is commented, then add it
			if((iContent.Key.empty())) {										// Is this a section?
				s += '[' + iContent.Section + ']';						// Add the section
			}
			else { s += iContent.Key + '=' + iContent.Value; }				// Or the Key value to the return srting

			i++;
			if (i != static_cast<int>(content.size())) { s += '\n'; }								// If this is not the last line, add a CrLf
		}
		return s;															// Return the contents
	}

	return "";
}

std::vector<std::string> CIniFile::GetSectionNames(const std::string& FileName)
{
	std::vector<std::string> data;													// Holds the return data
	std::vector<Record> content;													// Holds the current record													// Holds the current record

	if (Load(FileName, content))											// Make sure the file is loaded
	{
		for (const auto& iContent : content)								// Loop through the content
		{
			if(iContent.Key.empty()) {											// If there is no key value, then its a section
				data.push_back(iContent.Section);							// Add the section to the return data
			}
		}
	}

	return data;															// Return the data
}

std::vector<CIniFile::Record> CIniFile::GetSection(const std::string& SectionName, const std::string& FileName)
{
	std::vector<Record> data;													// Holds the return data
	std::vector<Record> content;													// Holds the current record													// Holds the current record

	if (Load(FileName, content))											// Make sure the file is loaded
	{
		for (const auto& iContent : content)								// Loop through the content
		{
			if((iContent.Section == SectionName) &&						// If this is the section name we want
				(!iContent.Key.empty()))										// but not the section name itself
				{ data.push_back(iContent); }									// Add the record to the return data
		}
	}
	
	return data;															// Return the data
}

bool CIniFile::RecordExists(std::string& KeyName, std::string& SectionName, const std::string& FileName)
{
	std::vector<Record> content;													// Holds the current record													// Holds the current record

	if (Load(FileName, content))											// Make sure the file is loaded
	{
		auto iter = std::find_if(content.begin(), 
				content.end(), 
				CIniFile::RecordSectionKeyIs(SectionName,KeyName));			// Locate the Section/Key

		if (iter == content.end()) { return false; }							// The Section/Key was not found
	}
	return true;															// The Section/Key was found
}

bool CIniFile::SectionExists(std::string& SectionName, const std::string& FileName)
{
	std::vector<Record> content;													// Holds the current record													// Holds the current record

	if (Load(FileName, content))											// Make sure the file is loaded
	{
		auto iter = std::find_if(content.begin(), 
				content.end(), 
				CIniFile::RecordSectionIs(SectionName));					// Locate the Section

		if (iter == content.end()) { return false; }							// The Section was not found
	}
	return true;															// The Section was found
}

std::vector<CIniFile::Record> CIniFile::GetRecord(std::string& KeyName, std::string& SectionName, const std::string& FileName)
{
	std::vector<Record> data;													// Holds the return data
	std::vector<Record> content;													// Holds the current record													// Holds the current record

	if (Load(FileName, content))											// Make sure the file is loaded
	{
		auto iter = std::find_if(content.begin(), 
				content.end(), 
				CIniFile::RecordSectionKeyIs(SectionName,KeyName));			// Locate the Record

		if (iter == content.end()) { return data; }								// The Record was not found

		data.push_back (*iter);												// The Record was found
	}
	return data;															// Return the Record
}

std::string CIniFile::GetValue(std::string& KeyName, std::string& SectionName, const std::string& FileName)
{
	if (!EERAModel::Utilities::fileExists(FileName)) { throw EERAModel::IO::IOException(FileName + ": File not found!"); }

	std::vector<Record> content = GetRecord(KeyName, SectionName, FileName);		// Get the Record

	if(!content.empty()) {													// Make sure there is a value to return
		return content[0].Value;											// And return the value
	}

	return "";																// No value was found
}


bool CIniFile::SetValue(std::string& KeyName, std::string& Value, std::string& SectionName, const std::string& FileName)
{
	std::vector<Record> content;													// Holds the current record													// Holds the current record

	if (Load(FileName, content))											// Make sure the file is loaded
	{
		if(!SectionExists(SectionName,FileName))							// If the Section doesn't exist
		{
			Record s = {"", ' ', SectionName, "", ""};							// Define a new section
			Record r = {"", ' ', SectionName, KeyName, Value};					// Define a new record
			content.push_back(s);											// Add the section
			content.push_back(r);											// Add the record
			return Save(FileName,content);									// Save
		}

		if(!RecordExists(KeyName,SectionName,FileName))						// If the Key doesn't exist
		{
			auto iter = std::find_if(content.begin(), 
				content.end(), 
				CIniFile::RecordSectionIs(SectionName));					// Locate the Section
			iter++;															// Advance just past the section
			Record r = {"",' ',SectionName,KeyName,Value};						// Define a new record
			content.insert(iter,r);											// Add the record
			return Save(FileName,content);									// Save
		}

		auto iter = std::find_if(content.begin(), 
				content.end(), 
				CIniFile::RecordSectionKeyIs(SectionName,KeyName));			// Locate the Record

		iter->Value = Value;												// Insert the correct value
		return Save(FileName,content);										// Save
	}

	return false;															// In the event the file does not load
}

bool CIniFile::RenameSection(std::string& OldSectionName, std::string& NewSectionName, const std::string& FileName)
{
	std::vector<Record> content;													// Holds the current record													// Holds the current record

	if (Load(FileName, content))											// Make sure the file is loaded
	{
		for(auto iter = content.begin(); 
			iter < content.end(); iter++)									// Loop through the records
		{
			if(iter->Section == OldSectionName)	{							// Is this the OldSectionName?
				iter->Section = NewSectionName;								// Now its the NewSectionName
			}
		}
		return Save(FileName,content);										// Save
	}

	return false;															// In the event the file does not load
}

bool CIniFile::CommentRecord(CommentChar& cc, std::string& KeyName, std::string& SectionName, const std::string& FileName)
{
	std::vector<Record> content;													// Holds the current record													// Holds the current record

	if (Load(FileName, content))											// Make sure the file is loaded
	{
		auto iter = std::find_if(content.begin(), 
				content.end(), 
				CIniFile::RecordSectionKeyIs(SectionName,KeyName));			// Locate the Section/Key

		if (iter == content.end()) { return false; }							// The Section/Key was not found
	
		iter->Commented = cc;										// Change the Comment value
		return Save(FileName,content);										// Save

	}
	return false;															// In the event the file does not load
}

bool CIniFile::UnCommentRecord(std::string& KeyName, std::string& SectionName, const std::string& FileName)
{
	std::vector<Record> content;													// Holds the current record													// Holds the current record

	if (Load(FileName, content))											// Make sure the file is loaded
	{
		auto iter = std::find_if(content.begin(), 
				content.end(), 
				CIniFile::RecordSectionKeyIs(SectionName,KeyName));			// Locate the Section/Key

		if (iter == content.end()) { return false; }							// The Section/Key was not found
	
		iter->Commented = ' ';												// Remove the Comment value
		return Save(FileName,content);										// Save

	}
	return false;															// In the event the file does not load
}

bool CIniFile::CommentSection(char& CommentChar, std::string& SectionName, const std::string& FileName)
{
	std::vector<Record> content;													// Holds the current record													// Holds the current record

	if (Load(FileName, content))											// Make sure the file is loaded
	{
		for(auto iter = content.begin(); iter < content.end(); iter++)
		{
			if(iter->Section == SectionName) {								// Is this the right section?
				iter->Commented = CommentChar;								// Change the comment value
			}
		}
		return Save(FileName,content);										// Save
	}

	return false;															// In the event the file does not load
}

bool CIniFile::UnCommentSection(std::string& SectionName, const std::string& FileName)
{
	std::vector<Record> content;													// Holds the current record													// Holds the current record

	if (Load(FileName, content))											// Make sure the file is loaded
	{
		for(auto iter = content.begin(); iter < content.end(); iter++)
		{
			if(iter->Section == SectionName) {								// Is this the right section?
				iter->Commented = ' ';										// Remove the comment value
			}
		}																	
		return Save(FileName,content);										// Save
	}

	return false;															// In the event the file does not load
}

bool CIniFile::DeleteRecord(std::string& KeyName, std::string& SectionName, const std::string& FileName)
{
	std::vector<Record> content;													// Holds the current record													// Holds the current record

	if (Load(FileName, content))											// Make sure the file is loaded
	{
		auto iter = std::find_if(content.begin(), 
				content.end(), 
				CIniFile::RecordSectionKeyIs(SectionName,KeyName));			// Locate the Section/Key

		if (iter == content.end()) { return false; }							// The Section/Key was not found
	
		content.erase(iter);												// Remove the Record
		return Save(FileName,content);										// Save

	}
	
	return false;															// In the event the file does not load
}

bool CIniFile::DeleteSection(std::string& SectionName, const std::string& FileName)
{
	std::vector<Record> content;													// Holds the current record													// Holds the current record

	if (Load(FileName, content))											// Make sure the file is loaded
	{
		for(int i=static_cast<int>(content.size()) - 1; i > -1; i--)								// Iterate backwards through the content
		{							
			if(content[i].Section == SectionName) {							// Is this related to the Section?
				content.erase (content.begin()+i);							// Then erase it
			}
		}

		return Save(FileName,content);										// Save
	}
	return false;															// In the event the file does not load
}

bool CIniFile::SetSectionComments(std::string& Comments, std::string& SectionName, const std::string& FileName)
{
	std::vector<Record> content;													// Holds the current record													// Holds the current record

	if (Load(FileName, content))											// Make sure the file is loaded
	{
		for(auto iter = content.begin(); iter < content.end(); iter++)									// Loop through the records
		{
			if((iter->Section == SectionName) &&							// Is this the Section?
				(iter->Key.empty()))											// And not a record
			{	
				if (Comments.size() >= 2)									// Is there a comment?
				{
					if (Comments.substr(Comments.size()-2) != "\n") {		// Does the string end in a newline?
						Comments += "\n";								// If not, add one
					}
				}
				iter->Comments = Comments;								// Set the comments
					
				return Save(FileName,content);							// Save
			}
		}
	}
	return false;															// In the event the file does not load
}

bool CIniFile::SetRecordComments(std::string& Comments, std::string& KeyName, std::string& SectionName, const std::string& FileName)
{
	std::vector<Record> content;													// Holds the current record													// Holds the current record

	if (Load(FileName, content))											// Make sure the file is loaded
	{
		auto iter = std::find_if(content.begin(), 
				content.end(), 
				CIniFile::RecordSectionKeyIs(SectionName,KeyName));			// Locate the Section/Key

		if (iter == content.end()) { return false; }							// The Section/Key was not found
	
		if (Comments.size() >= 2)											// Is there a comment?
		{
			if (Comments.substr(Comments.size()-2) != "\n")	{				// Does the string end in a newline?
				Comments += "\n";											// If not, add one
			}
		}
		iter->Comments = Comments;											// Set the comments
		return Save(FileName,content);										// Save

	}
	
	return false;															// In the event the file does not load
}

std::vector<CIniFile::Record> CIniFile::GetSections(const std::string& FileName)
{
	std::vector<Record> data;													// Holds the return data
	std::vector<Record> content;													// Holds the current record

	if (Load(FileName, content))											// Make sure the file is loaded
	{
		for (const auto& iContent : content)								// Loop through the content
		{
			if(iContent.Key.empty()) {										// If this is a section 
				data.push_back(iContent);									// Add the record to the return data
			}
		}
	}
	
	return data;															// Return the data
}

bool CIniFile::AddSection(std::string& SectionName, const std::string& FileName)
{
	std::vector<Record> content;													// Holds the current record

	if (Load(FileName, content))											// Make sure the file is loaded
	{
		Record s = {"", ' ', SectionName, "", ""};								// Define a new section
		content.push_back(s);												// Add the section
		return Save(FileName,content);										// Save
	}

	return false;															// The file did not open
}

bool CIniFile::Create(const std::string& FileName)
{
	std::vector<Record> content;													// Create empty content
	return Save(FileName,content);											// Save
}

