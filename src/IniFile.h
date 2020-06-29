#pragma once
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <functional>


class CIniFile
{
public:
	struct Record
	{
		std::string Comments;
		char Commented;
		std::string Section;
		std::string Key;
		std::string Value;
	};

	enum CommentChar
	{
		Pound = '#',
		SemiColon = ';'
	};

	CIniFile();

	static bool AddSection(std::string& SectionName, const std::string& FileName);
	static bool CommentRecord(CommentChar& cc, std::string& KeyName, std::string& SectionName, const std::string& FileName);
	static bool CommentSection(char& CommentChar, std::string& SectionName, const std::string& FileName);
	static std::string Content(const std::string& FileName);
	static bool Create(const std::string& FileName);
	static bool DeleteRecord(std::string& KeyName, std::string& SectionName, const std::string& FileName);
	static bool DeleteSection(std::string& SectionName, const std::string& FileName);
	static std::vector<Record> GetRecord(std::string& KeyName, std::string& SectionName, const std::string& FileName);
	static std::vector<Record> GetSection(const std::string& SectionName, const std::string& FileName);
	static std::vector<std::string> GetSectionNames(const std::string& FileName);
	static std::string GetValue(std::string& KeyName, std::string& SectionName, const std::string& FileName);
	static bool RecordExists(std::string& KeyName, std::string& SectionName, const std::string& FileName);
	static bool RenameSection(std::string& OldSectionName, std::string& NewSectionName, const std::string& FileName);
	static bool SectionExists(std::string& SectionName, const std::string& FileName);
	static bool SetRecordComments(std::string& Comments, std::string& KeyName, std::string& SectionName, const std::string& FileName);
	static bool SetSectionComments(std::string& Comments, std::string& SectionName, const std::string& FileName);
	static bool SetValue(std::string& KeyName, std::string& Value, std::string& SectionName, const std::string& FileName);
	static bool UnCommentRecord(std::string& KeyName, std::string& SectionName, const std::string& FileName);
	static bool UnCommentSection(std::string& SectionName, const std::string& FileName);

private:
	static std::vector<Record> GetSections(const std::string& FileName);
	static void Trim(std::string& str, const std::string& CharsToTrim = " \t\n\r", int TrimDir = 0);
	static bool Load(const std::string& FileName, std::vector<Record>& content);	
	static bool Save(const std::string& FileName, const std::vector<Record>& content);

	struct RecordSectionIs : std::unary_function<Record, bool>
	{
		std::string section_;

		explicit RecordSectionIs(std::string& section): section_(section){}

		bool operator()( const Record& rec ) const
		{
			return rec.Section == section_;
		}
	};

	struct RecordSectionKeyIs : std::unary_function<Record, bool>
	{
		std::string section_;
		std::string key_;

		explicit RecordSectionKeyIs(std::string& section, std::string& key): section_(section),key_(key){}

		bool operator()( const Record& rec ) const
		{
			return ((rec.Section == section_)&&(rec.Key == key_));
		}
	};

	struct AscendingSectionSort
	{
		bool operator()(Record& Start, Record& End)
		{
			return Start.Section < End.Section;
		}
	};

	struct DescendingSectionSort
	{
		bool operator()(Record& Start, Record& End)
		{
			return Start.Section > End.Section;
		}
	};

	struct AscendingRecordSort
	{
		bool operator()(Record& Start, Record& End)
		{
			return Start.Key < End.Key;
		}
	};

	struct DescendingRecordSort
	{
		bool operator()(Record& Start, Record& End)
		{
			return Start.Key > End.Key;
		}
	};
};
