/*******************************************************************
this is adapted from the open source project: simple access layer
https://github.com/simple-access-layer/source
European Union Public License 1.2
https://github.com/simple-access-layer/source/blob/master/LICENCE.txt
by Qingfeng Xia UKAEA, 2020
*******************************************************************/
#pragma once

#include <string>
#include <memory>

#include "json.hpp"

#include "Poco/Net/HTTPRequest.h"
#include "Poco/Net/HTTPResponse.h"
#include "Poco/Net/HTTPSClientSession.h"
#include "Poco/Net/HTTPClientSession.h"
#include "Poco/StreamCopier.h"
#include "Poco/URI.h"
#include "Poco/UTF8String.h"

//#include "Poco/JSON/Object.h"
//#include "Poco/JSON/Parser.h"
//#include "Poco/SharedPtr.h"

/// utility class to hold static methods that work in header-only mode
class Restful
{
public:
    using JSON = nlohmann::json;

    /// there is no auth, so it will only work within intranet
    static JSON get_json(const std::string url)
    {
        using namespace Poco;
        using namespace Poco::Net;

        HTTPResponse response;
        HTTPClientSession *session;
        std::stringstream ss;
        Poco::URI uri{url};

        if (url.find("https://") != std::string::npos)
        {
            // https session, ignoring invalid certificates
            session = new HTTPSClientSession(uri.getHost(), uri.getPort(),
                                             new Context(Context::CLIENT_USE, "", Context::VERIFY_NONE));
        }
        else
        {
            // http session
            session = new HTTPClientSession(uri.getHost(), uri.getPort());
        }
        HTTPRequest request(HTTPRequest::HTTP_GET, url, HTTPMessage::HTTP_1_1);
        session->sendRequest(request);
        size_t bs = StreamCopier::copyStream64(session->receiveResponse(response), ss);
        // cout << "res.getStatus() = " << response.getStatus() << std::endl;

        // Poco's JSON api is hard to use IMO
        //using namespace Poco::JSON;
        //Poco::JSON::Parser parser;
        //auto j = parser.parse(json).extract<Poco::JSON::Object::Ptr>(); // decode json

        return JSON::parse(ss.str());
    }
};