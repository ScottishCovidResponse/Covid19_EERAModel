// Based on an example from https://github.com/andrew-hardin/cmake-git-version-tracking

// MIT License

// Copyright (c) 2020 Andrew Hardin

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#pragma once

#include <string>

class GitMetadata {
public:
  // Is the metadata populated? We may not have metadata if
  // there wasn't a .git directory (e.g. downloaded source
  // code without revision history).
  static bool Populated();

  // Were there any uncommitted changes that won't be reflected
  // in the CommitID?
  static bool AnyUncommittedChanges();

  // The commit author's name.
  static std::string AuthorName();
  // The commit author's email.
  static std::string AuthorEmail();
  // The commit SHA1.
  static std::string CommitSHA1();
  // The ISO8601 commit date.
  static std::string CommitDate();
  // The commit subject.
  static std::string CommitSubject();
  // The commit body.
  static std::string CommitBody();
  // The commit tag
  static std::string Tag();
};
